The Standard Gate Set
=====================

Below is a summary of the key gates used in tweedledum

.. |id_matrix| replace:: :math:`\pmatrix{1&0 \\ 0&1}`
.. |h_matrix|  replace:: :math:`\frac{1}{\sqrt{2}}\pmatrix{1&1 \\ 1&-1}`

.. |x_matrix| replace:: :math:`\pmatrix{0&1 \\ 1&0}`
.. |z_matrix| replace:: :math:`\pmatrix{1&0 \\ 0&-1}`

.. |t_matrix| replace:: :math:`\pmatrix{1&0 \\ 0&e^{\mathrm{i}\frac\pi4}}`
.. |s_matrix| replace:: :math:`\pmatrix{1&0 \\ 0&\mathrm{i}}`

.. |rz_matrix| replace:: :math:`\pmatrix{e^{-\mathrm{i}\theta}&0 \\ 0&e^{\mathrm{i}\theta}}`
.. |rx_matrix| replace:: :math:`\pmatrix{\cos\frac\theta2 & -\mathrm{i}\sin\frac\theta2 \\ -\mathrm{i}\sin\frac\theta2 & \cos\frac\theta2}`

.. |cx_matrix| replace:: :math:`\pmatrix{1&0&0&0 \\ 0&1&0&0 \\ 0&0&0&1 \\ 0&0&1&0}`
.. |cz_matrix| replace:: :math:`\pmatrix{1&0&0&0 \\ 0&1&0&0 \\ 0&0&1&0 \\ 0&0&0&-1}`

+--------------------------------+--------+---------------------------+---------------+
| Name(s)                        | Symbol | tweedledum symbol         |  Matrix       |
+=================+==============+========+===========================+===============+
| Identity                       | I      | ``gate_set::identity``    | |id_matrix|   |
+--------------------------------+--------+---------------------------+---------------+
| Hadamard                       | H      | ``gate_set::hadamard``    | |h_matrix|    |
+--------------------------------+--------+---------------------------+---------------+
| .. centered:: **Arbitrary rotations**                                               |
+--------------------------------+--------+---------------------------+---------------+
| X Rotation                     | Rx     | ``gate_set::rotation_x``  | |rx_matrix|   |
+--------------------------------+--------+---------------------------+---------------+
| Y Rotation                     | Ry     | ``gate_set::rotation_y``  |               |
+--------------------------------+--------+---------------------------+---------------+
| Z Rotation                     | Rz     | ``gate_set::rotation_z``  | |rz_matrix|   |
+--------------------------------+--------+---------------------------+---------------+
| .. centered:: **Named Rotations**                                                   |
+--------------------------------+--------+---------------------------+---------------+
| Pauli X, NOT                   | X      | ``gate_set::pauli_x``     | |x_matrix|    |
+--------------------------------+--------+---------------------------+---------------+
| T                              | T      | ``gate_set::t``           | |t_matrix|    |
+--------------------------------+--------+---------------------------+---------------+
| T dagger                       | T†     | ``gate_set::t_dagger``    |               |
+--------------------------------+--------+---------------------------+---------------+
| Phase                          | S      | ``gate_set::phase``       | |s_matrix|    |
+--------------------------------+--------+---------------------------+---------------+
| Phase dagger                   | S†     | ``gate_set::phase_dagger``|               |
+--------------------------------+--------+---------------------------+---------------+
| Pauli Z, Phase flip            | Z      | ``gate_set::pauli_z``     | |z_matrix|    |
+--------------------------------+--------+---------------------------+---------------+
| .. centered:: **Controlled gates**                                                  |
+--------------------------------+--------+---------------------------+---------------+
| Control NOT                    | CNOT   | ``gate_set::cx``          | |cx_matrix|   |
+--------------------------------+--------+---------------------------+---------------+
| Control Z                      | CZ     | ``gate_set::cz``          | |cz_matrix|   |
+--------------------------------+--------+---------------------------+---------------+
| Multiple Control NOT, Toffoli  |        | ``gate_set::mcx``         |               |
+--------------------------------+--------+---------------------------+---------------+
| Multiple Control Z             |        | ``gate_set::mcz``         |               |
+--------------------------------+--------+---------------------------+---------------+
