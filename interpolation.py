###lagrange interpolation

def multiplypoly(poly1,poly2):
    result=[0]*(len(poly1)+len(poly2)-1)
    for i in range(len(poly1)):
        for j in range(len(poly2)):
            result[i+j]+=poly1[i]*poly2[j]
    return result

def computeLk(points,k):
    n=len(points)
    x_points=[p[0] for p in points]
    
    num=[1.0]
    den=1
    for j in range(n):
            if j!=k:
                num=multiplypoly(num,[-x_points[j],1])
                den*=(x_points[k]-x_points[j])
    return [c/den for c in num]


def computePL(points):
    n=len(points)
    x_points=[p[0] for p in points]
    y_points=[p[1] for p in points]
    p=[0.0]*n
    for k in range(n):
        Lk=computeLk(points,k)
        for i in range(len(Lk)):
            p[i]+=y_points[k]*Lk[i]

    return p
    



def all_Lk(points):
    n=len(points)
    
    for k in range(n):
        Lk=computeLk(points,k)

        print(f"L_{k}(x)=", end="")
        style=" ¹²³⁴⁵⁶⁷⁸⁹"
        for d in range(len(Lk)-1,-1,-1):
            if round(Lk[d], 2) == 0:
                continue
            if (d==0) :
                print(f"{Lk[d]:.2f}",end=" ")
            else:
                print(f"{Lk[d]:.2f}x{style[d]}",end=" + ")
        print()


   


def lagrange_interpolation(points):
   all_Lk(points)
   print()
   computePL(points)


#### newton polynomial



def compute_Nk(points,k):
    x_points=[p[0] for p in points]
    
    num=[1.0]
    for i in range(k):
        num=multiplypoly(num,[-x_points[i],1.0])
    return num
     


def coefficients(points):
    x_points = [p[0] for p in points]
    y_points = [p[1] for p in points]
    n = len(points)
    
    if n == 0:
        return []
    
    # Copie des ordonnées 
    coeffs = [0.0] * n
    for i in range(n):
        coeffs[i] = y_points[i]
    
   
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (x_points[i] - x_points[i - j])
            
    return coeffs




def compute_pN(points):
    n=len(points)
    coeffs=coefficients(points)
    result=[0.0]*(n)
    for i in range(n):
        Nk=compute_Nk(points,i)
        for j in range(len(Nk)):
            result[j]+=coeffs[i]*Nk[j]

    return result
#### hermite interpolation


def compute_Ĥk(points,k):

    x_points=[p[0] for p in points]
    # 1. Calcul du polynôme de Lagrange L_k
    Lk = computeLk(points, k)
    
    # 2. Calcul de L_k^2
    Lk_squared = multiplypoly(Lk, Lk)
    
    # 3. Multiplication par (x - x_k)
    factor = [-x_points[k], 1]
    return multiplypoly(factor, Lk_squared)

def derive_lK(poly):
    derived = []
    
    for i in range(1, len(poly)):
        derived = derived + [i * poly[i]]
        
    return derived

def compute_Hk(points,k):
    x_points=[p[0] for p in points]
        #  Calcul du polynôme de Lagrange L_k
    Lk = computeLk(points, k)
        
        # Calcul de L_k^2
    Lk_squared = multiplypoly(Lk, Lk)

        # calcul de la derivée de L_k
    Lk_derivative = derive_lK(Lk)

    def evaluate_derivative_at_xk(Lk_derivative,x_k):
        
        result = 0
        for i in range(len(Lk_derivative)):
            result += Lk_derivative[i] * (x_k ** i)
        return result

    
    factor = [1 + 2 * x_points[k] * evaluate_derivative_at_xk(Lk_derivative, x_points[k]), -2 * evaluate_derivative_at_xk(Lk_derivative, x_points[k])]

    return multiplypoly(factor, Lk_squared)

def compute_pH(points):
    n=len(points)
    x_points=[p[0] for p in points]
    y_points=[p[1] for p in points]
    derivatives=[p[2] for p in points]
    p=[0]*(2*n)

    for k in range(n):
        Hk = compute_Hk(points, k)
        Ĥk = compute_Ĥk(points, k)
        for i in range(len(Hk)):
            p[i] += y_points[k] * Hk[i]
            p[i] += derivatives[k] * Ĥk[i]

    return p
def evaluer_polynome(coeffs, x):

    valeur = 0.0
    for i, c in enumerate(coeffs):
        valeur += c * (x ** i)
    return valeur
