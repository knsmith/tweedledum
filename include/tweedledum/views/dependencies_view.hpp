/*-------------------------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Bruno Schmitt
*------------------------------------------------------------------------------------------------*/
#pragma once

#include "../traits.hpp"
#include "../utils/node_map.hpp"
#include "immutable_view.hpp"

#include <array>
#include <vector>

namespace tweedledum {

/*! \brief This view computes real dependencies between pairs of consective gates.
 *
 *  In a gate graph network (`gg_network`) connected pair of nodes have a parent-child dependency
 *  relation. However, if a pair o gates commutes, the dependency is false and these gate do not
 *  depend on each other, i.e., they can apper on either order. 
 * 
 *  The view modifies the behavior of `foreach_child` and implements get_dependencies method. 
 *  
 * **Required gate functions:**
 * - `is_dependent`
 * - `foreach_control`
 * - `foreach_target`
 * - `qubit_slot`
 *
 * **Required network functions:**
 * - `foreach_cgate`
 * - `foreach_cinput`
 * - `foreach_coutput`
 * - `get_node`
 * - `get_output`
 * 
 */
template<typename Network>
class dependencies_view : public immutable_view<Network> {
public:
	using gate_type = typename Network::gate_type;
	using node_type = typename Network::node_type;
	using node_ptr_type = typename Network::node_ptr_type;
	using storage_type = typename Network::storage_type;

	/*! \brief Default constructor.
	 *
	 * Constructs depth view on a network.
	 */
	explicit dependencies_view(Network& network)
	    : immutable_view<Network>(network)
	    , dependencies_(network)
	{
		compute_dependencies();
	}

#pragma region Node properties
	/*! \brief Returns the node dependencies. */
	auto const& get_dependencies(node_type const& node) const
	{
		return dependencies_[node];
	}

	/*! \brief Returns the node dependencies w.r.t a qubit. */
	auto const& get_dependencies(node_type const& node, qubit_id qid) const
	{
		return dependencies_[node][node.gate.qubit_slot(qid)];
	}
#pragma endregion

#pragma region Const iterators
	/*! \brief Calls ``fn`` on all **dependent** children of a node.
	 *
	 * The paramater ``fn`` is any callable that must have one of the following two signatures:
	 * - ``void(node_ptr const&)``
	 * - ``void(node_ptr const&, node_ptr_type)``
	 * 
	 * If ``fn`` has two parameters, the second parameter is an index starting
	 * from 0 and incremented in every iteration.
	 */
	template<typename Fn>
	void foreach_child(node_type const& node, Fn&& fn) const
	{
		// clang-format off
		static_assert(is_callable_without_index_v< Fn, node_ptr_type, void>
		              || is_callable_with_index_v<Fn, node_ptr_type, void>);
		// clang-format on

		auto& dependent_nodes = dependencies_[node];
		for (auto i = 0u; i < dependent_nodes.size(); ++i) {
			for (auto child_index : dependent_nodes[i]) {
				if constexpr (is_callable_without_index_v<Fn, node_ptr_type, void>) {
					fn(child_index);
				} else if constexpr (is_callable_with_index_v<Fn, node_ptr_type, void>) {
					fn(child_index, i);
				}
			}
		}
	}

	/*! \brief Calls ``fn`` on all **dependent** children w.r.t a qubit of a node.
	 *
	 * The paramater ``fn`` is any callable that must have one of the following signature:
	 * - ``void(node_ptr const&)``
	 */
	template<typename Fn>
	void foreach_child(node_type const& node, qubit_id qid, Fn&& fn) const
	{
		auto const& node_dependency = get_dependencies(node, node.gate.qubit_slot(qid));
		for (auto dependency_index : node_dependency) {
			fn(dependency_index);
		}
	}
#pragma endregion

private:
#pragma region Non-const iterators
	template<typename Fn>
	void foreach_child(node_type const& node, qubit_id qid, Fn&& fn)
	{
		auto& node_dependency = dependencies_[node].at(node.gate.qubit_slot(qid));
		for (auto dependency_index : node_dependency) {
			fn(dependency_index);
		}
	}
#pragma endregion

	void connect_node(qubit_id qid, node_type const& node)
	{
		const auto& output = this->get_output(qid);
		auto& output_dependency = dependencies_[output].at(output.gate.qubit_slot(qid));
		auto& node_dependency = dependencies_[node].at(node.gate.qubit_slot(qid));
		// Previous node
		auto prev_node_index = output_dependency.back();
		auto& prev_node = this->get_node(prev_node_index);

		if (node.gate.is_dependent(prev_node.gate)) {
			for (auto index : output_dependency) {
				node_dependency.emplace_back(index);
			}
			output_dependency.clear();
			output_dependency.push_back(this->node_to_index(node));
			return;
		}
		output_dependency.push_back(this->node_to_index(node));
		foreach_child(prev_node, qid, [&](auto arc) {
			node_dependency.emplace_back(arc);
		});
	}

	void compute_dependencies()
	{
		// Start by connecting input and ouput nodes of the same qubit
		this->foreach_cinput([&](auto const& node, auto node_index) {
			const auto& gate = node.gate;
			gate.foreach_target([&](auto qid) {
				const auto& output = this->get_output(qid);
				auto& dependency = dependencies_[output].at(gate.qubit_slot(qid));
				dependency.push_back(node_index);
			});
		});

		this->foreach_cgate([&](auto const& node) {
			node.gate.foreach_control([&](auto qid) { connect_node(qid, node); });
			node.gate.foreach_target([&](auto qid) { connect_node(qid, node); });
		});
	}

private:
	using dependecy_array = std::array<std::vector<node_ptr_type>, gate_type::max_num_qubits>;
	node_map<dependecy_array, Network> dependencies_;
};

} // namespace tweedledum
