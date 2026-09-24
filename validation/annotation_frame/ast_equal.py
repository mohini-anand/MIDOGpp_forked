"""
    Prove a docstring-only edit is inert: two Python files are equal once every docstring is stripped from their ASTs.
    Usage: ast_equal.py A B
"""

import ast
import sys


def stripped(path):
    """
        The AST dump of a file with every module, class and function docstring removed.

        path (str): the Python file.

        Returns str: the dump.
    """
    tree = ast.parse(open(path).read())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:]
    return ast.dump(tree)


if __name__ == "__main__":
    a, b = sys.argv[1], sys.argv[2]
    same = stripped(a) == stripped(b)
    print(f"docstring-stripped ASTs equal: {same} ({a} vs {b})")
    sys.exit(0 if same else 1)
