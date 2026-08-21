import math
import numpy as np


def calculate_reverse_triangle_matrix(A):
    """Evaluates the matrix formulation of the reverse triangle equation:

    ((P_g * M * P_g^T) * (Q_g * L))_{g=0}^T

    Parameters:
        A (list or array): Input array of pyramid items before differences.

    Returns:
        dict: A dictionary containing:
            - 'results': The outer list of matrix outputs matching the iterative logic
            - 'T': The calculated value of T
            - 'terms_needed': The amount of terms needed to calculate
            (tetrahedral number)
    """
    T = len(A) - 1

    if T < 0:
        raise ValueError("Input array A must contain at least one item.")

    # 1. Construct Matrix M of size (T + 1) x (T + 1)
    # M_{k,n} = (-1)^(n+k) * binom(k, n) if 0 <= n <= k <= T else 0
    M = np.zeros((T + 1, T + 1), dtype=int)
    for k in range(T + 1):
        for n in range(k + 1):
            M[k, n] = ((-1) ** (n + k)) * math.comb(k, n)

    # 2. Construct Column Vector L of size (T + 1) x 1
    # L_{J, 0} = A[J]
    L = np.array(A, dtype=int).reshape((T + 1, 1))

    outer_results = []

    # Loop for g from 0 to T
    for g in range(T + 1):
        dim = T + 1 - g  # Size: (T + 1 - g)

        # 3. Construct Projection Matrix P_g of size (T + 1 - g) x (T + 1)
        # P_{g; i, j} = 1 if 0 <= i = j <= T - g else 0
        P_g = np.zeros((dim, T + 1), dtype=int)
        for i in range(dim):
            P_g[i, i] = 1

        P_g_T = P_g.T  # Transpose of P_g, size (T + 1) x (T + 1 - g)

        # 4. Construct Shift/Selection Matrix Q_g of size (T + 1 - g) x (T + 1)
        # (Q_g)_{i, j} = 1 if j = i + g else 0
        Q_g = np.zeros((dim, T + 1), dtype=int)
        for i in range(dim):
            Q_g[i, i + g] = 1

        # Matrix multiplications for current index g:
        # Step A: Reduced M matrix -> (P_g @ M @ P_g_T) of size (dim x dim)
        M_reduced = P_g @ M @ P_g_T

        # Step B: Shifted/Sliced L vector -> (Q_g @ L) of size (dim x 1)
        L_shifted = Q_g @ L

        # Step C: Multiply reduced M with shifted L -> size (dim x 1)
        res_vector = M_reduced @ L_shifted

        # Flatten column vector into list format for output matching
        outer_results.append(res_vector.flatten().tolist())

    terms_needed = math.comb(T + 3, 3)

    return {"results": outer_results, "T": T, "terms_needed": terms_needed}


# --- Example Usage ---
if __name__ == "__main__":
    A_input = [10, 20, 35, 60]

    output = calculate_reverse_triangle_matrix(A_input)

    print(f"Input Array A: {A_input}")
    print(f"T (len(A) - 1): {output['T']}")
    print(
        f"Amount of terms needed (Tetrahedral Number): {output['terms_needed']}"
    )
    print("\nMatrix Evaluated Output:")
    for row in output["results"]:
        print(row)
