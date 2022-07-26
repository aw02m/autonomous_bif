import numpy as np
import sympy as sp
import os

# # Rossler
# xdim = 3
# pdim = 3

# def f(x, p):
#     return sp.Matrix([
#         -x[1] - x[2],
#         x[0] + p[0] * x[1],
#         p[1] * x[0] - p[2] * x[2] + x[0] * x[2]
#     ])

# coupled FHN
# xdim: dimension of the system
# pdim: dimension of the parameter
xdim = 4
pdim = 5


def h(x):
    return (1 + sp.tanh(x)) / 2


def func(x, p):  # define the function here
    return sp.Matrix([
        p[4] * (x[0] - p[2] * x[0] * x[0] * x[0] - x[2] + p[3] * h(x[1])),
        p[4] * (x[1] - p[2] * x[1] * x[1] * x[1] - x[3] - p[3] * h(x[0])),
        x[0] + p[0] - p[1] * x[2],
        x[1] + p[0] - p[1] * x[3]
    ])

# DO NOT EDIT BELOW


# output path where the cpp function file is wrote
pp_output_path = '../pp/cmake-tree/src/sys_func.cpp'
bif_output_path = '../bif/cmake-tree/src/sys_func.cpp'

header = """\
#include "dynamical_system.hpp"

void dynamical_system::sys_func(const Eigen::VectorXd &x, const double /*t*/) {
"""

sym_x = sp.MatrixSymbol('x', xdim, 1)
sym_p = sp.MatrixSymbol('p', pdim, 1)

n = xdim

code = ""

f = func(sym_x, sym_p)
for i in range(n):
    idx = i
    code += sp.ccode(f[idx], assign_to=('f('+str(i)+')'), standard='C89')
    code += "\n"
code += "\n"

dfdx = sp.derive_by_array([f[i] for i in range(xdim)], [
                          sym_x[i] for i in range(xdim)]).transpose()
for i in range(n):
    for j in range(n):
        code += sp.ccode(dfdx[i][j], assign_to=(
            'dfdx('+str(i)+', '+str(j)+')'), standard='C89')
        code += "\n"
code += "\n"

code += "switch (var_param) {\n"
for idx_param in range(pdim):
    code += "case "+str(idx_param)+":\n"
    dfdlambda = sp.diff(f, sym_p[idx_param])
    for i in range(n):
        idx = i
        code += sp.ccode(dfdlambda[idx], assign_to=(
            'dfdlambda('+str(i)+')'), standard='C89')
        code += "\n"
    code += "break;\n"
code += "}\n"
code += "\n"

dfdxdx = [sp.zeros(xdim, xdim) for j in range(xdim)]
for i in range(xdim):
    dfdxdx[i] = sp.diff(dfdx, sym_x[i])
for i in range(n):
    for j in range(n):
        for k in range(n):
            code += sp.ccode(dfdxdx[i][j][k], assign_to=(
                'dfdxdx['+str(i)+']('+str(j)+', '+str(k)+')'), standard='C89')
            code += "\n"
code += "\n"

code += "switch (var_param) {\n"
for idx_param in range(pdim):
    code += "case "+str(idx_param)+":\n"
    dfdxdlambda = sp.diff(dfdx, sym_p[idx_param])
    for i in range(n):
        for j in range(n):
            code += sp.ccode(dfdxdlambda[i][j], assign_to=(
                'dfdxdlambda('+str(i)+', '+str(j)+')'), standard='C89')
            code += "\n"
    code += "break;\n"
code += "}\n\n"

footer = """\
}"""

file = open(bif_output_path, 'w')
file.write(header + code + footer)

file.close()

os.system('clang-format -i -style=file ' + bif_output_path)

# pp part

header = """\
#include "dynamical_system.hpp"

Eigen::VectorXd dynamical_system::func([[maybe_unused]] double t,
                                       const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(xdim);

"""

code = ""

for i in range(n):
    idx = i
    code += sp.ccode(f[idx], assign_to=('ret('+str(i)+')'), standard='C89')
    code += "\n"
code += "\n"

footer = """\
  return ret;
}"""

file = open(pp_output_path, 'w')
file.write(header + code + footer)

file.close()

os.system('clang-format -i -style=file ' + pp_output_path)
