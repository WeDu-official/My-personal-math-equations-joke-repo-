# by weduofficial/weduxox
import math
def calculate_reverse_triangle(A):
    """
    Evaluates the nested equation:
    (( sum_{n=0}^K (-1)^(n+K) * binom(K, n) * A[n+g] )_{K=0}^{T-g})_{g=0}^{T}
    
    Parameters:
        A (list or array): The array of items [x, y, ...]
        
    Returns:
        dict: A dictionary containing:
            - 'results': The nested evaluated output structure
            - 'T': The calculated value of T
            - 'terms_needed': The amount of terms needed to calculate (tetrahedral number)
    """
    T = len(A) - 1
    
    if T < 0:
        raise ValueError("Input array A must contain at least one item.")

    # Calculate tetrahedral number: binom(T + 3, 3) = (T + 1)(T + 2)(T + 3) / 6
    terms_needed = math.comb(T + 3, 3)
    
    outer_results = []
    
    # Outer loop: g from 0 to T
    for g in range(T + 1):
        inner_results = []
        
        # Middle loop: K from 0 to T - g
        for K in range(T - g + 1):
            
            # Inner summation: n from 0 to K
            sum_val = 0
            for n in range(K + 1):
                sign = (-1) ** (n + K)
                binom_coeff = math.comb(K, n)
                term = sign * binom_coeff * A[n + g]
                sum_val += term
                
            inner_results.append(sum_val)
            
        outer_results.append(inner_results)
        
    return {
        "results": outer_results,
        "T": T,
        "terms_needed": terms_needed
    }

# --- Example Usage ---
if __name__ == "__main__":
    # Input array A (e.g., [x, y, z, w])
    A_input = [10, 20, 35, 60]
    
    output = calculate_reverse_triangle(A_input)
    
    print(f"Input Array A: {A_input}")
    print(f"T (len(A) - 1): {output['T']}")
    print(f"Amount of terms needed (Tetrahedral Number): {output['terms_needed']}")
    print("Nested Result Matrix:")
    for row in output["results"]:
        print(row)
