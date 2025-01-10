
import sympy
import numpy as np


# def get
from matplotlib import colormaps
from matplotlib import pyplot as plt
# def find_rest_vv


import sympy



"""
Everything here is adapted from the code in this notebook:
2024_11_27_sijie_system.ipynb
which repeated the previous caculations I did in terms of the energy paraameteers used by Sijie in his paper,
and by most people working with vertex models.
"""
class RegHexagonalModel:
    @classmethod
    def create_area_symexp(
        cls,
        s_alpha,
        s_beta,
    ):
        s_A = (sympy.sqrt(3)*3/2)*s_alpha*s_beta
        return s_A

    @classmethod
    def create_perim_symexp(
        cls,
        s_alpha,
        s_beta,
    ):
        s_P = 2*s_alpha + 2*sympy.sqrt(s_alpha**2 + 3 * s_beta**2)
        return s_P
    
    @classmethod
    def area_energy_symexp(
        cls,
        s_alpha,
        s_beta,
        s_Area_rest,
        s_K,
    ):
        s_Area = cls.create_area_symexp(s_alpha, s_beta)
        return (s_K / 2) * (s_Area - s_Area_rest)**2

        
    @classmethod
    def perimeter_energy_symexp(
        cls,
        s_alpha,
        s_beta,
        s_P_r,
        s_Gam,
    ):
        s_Perimeter = cls.create_perim_symexp(s_alpha, s_beta)
        return (s_Gam/2) * (s_Perimeter - s_P_r)**2
        
    
    @classmethod
    def create_energy_symexp(
        cls,
        s_alpha,
        s_beta,
        s_Area_rest,
        s_Perimeter_rest,
        s_K,
        s_Gam,
    ):
        s_E_area = cls.area_energy_symexp(s_alpha, s_beta, s_Area_rest, s_K)
        s_E_perim = cls.perimeter_energy_symexp(s_alpha, s_beta, s_Perimeter_rest, s_Gam)
        
        
        
        s_E = s_E_area + s_E_perim
        
        return s_E
        
    
    def __init__(self):
        alpha, beta = sympy.symbols("alpha beta", real=True, positive=True)
        A_r, P_r = sympy.symbols("A_r P_r", real=True, positive=True)
        K, Gam = sympy.symbols("K Gamma", real=True)
        mu = sympy.symbols("mu", real=True, positive=True)
        


        Gam_til, K_til, alpha_til =  sympy.symbols(
                                                        r"\tilde{\Gamma} \tilde{K} \tilde{\alpha}",
                                                        real=True,
                                                        positive=True)
                                                        
        beta_til, P_r_til, A_r_til =  sympy.symbols(
                                                        r"\tilde{\beta}  \tilde{P_r} \tilde{A_r}",
                                                        real=True,
                                                        positive=True)
        Y = sympy.symbols("Y", real=True, positive=True)
        xi_til =  sympy.symbols(r"\tilde{\xi}", real=True, positive=True)
        
        self.alpha = alpha
        self.beta = beta
        self.A_r = A_r
        self.P_r = P_r
        self.K = K
        self.Gam = Gam
        self.mu = mu
        self.Gam_til = Gam_til
        self.K_til = K_til
        self.alpha_til = alpha_til
        self.beta_til = beta_til
        self.P_r_til = P_r_til
        self.A_r_til = A_r_til
        self.Y = Y
        self.xi_til = xi_til


        self.energy_basic = self.create_energy_symexp(
            alpha,beta, A_r, P_r, K, Gam,
        )
        
                                                        
        self.tilde_mu_substitutions = [
            (alpha, alpha_til *mu),
            (beta, beta_til*mu),
            (P_r, P_r_til*mu),
            (Gam, Gam_til/mu**2),
            (A_r,  A_r_til*mu**2),
            (K, K_til / mu**4),
        ]
        
        
        self.E_mu_subbed = self.energy_basic.subs(
            self.tilde_mu_substitutions
        ).simplify().subs([
            (A_r_til, 1)
        ])
        
        
        self.E_til = (1/K_til) * self.E_mu_subbed.subs([
            ( Gam_til, Y * K_til ),
        ]).simplify()
        
        # Energy, if we scale vertically and horizontally by the same amount
        self.E_til_iso = self.E_til.subs([
            (alpha_til, xi_til),
            (beta_til, xi_til)
        ])
        

        self.dE_dxi_symexp = self.E_til_iso.diff(xi_til).expand()
        
        
    
        self.dE_dxi_polynomial = sympy.Poly(self.dE_dxi_symexp, xi_til)

        assert self.dE_dxi_polynomial.degree() == 3, "should be a cubic polynomial"
    
        # d2E_d2xi_polynomial
        self.d2E_dxi2_sym = self.dE_dxi_symexp.diff(xi_til).expand()
        self.lambified_d2e_dxi2 = sympy.lambdify((xi_til, Y,P_r_til), self.d2E_dxi2_sym)

        self.dE_dxi_coeff_functions = []
        for coeff_sym in self.dE_dxi_polynomial.all_coeffs():
            self.dE_dxi_coeff_functions.append(
                sympy.lambdify((Y,P_r_til), coeff_sym)
            )
        
        self.get_gradient, self.get_hessian = self.generate_lambdified_grad_and_hessian_function()
        
        
            
        
    def find_xi_til_equilibrium_numeric(self, Y_value, P_r_til_value):
        cubic_coeffs = [coeff_f(Y_value, P_r_til_value) for coeff_f in self.dE_dxi_coeff_functions]
        np_roots = np.roots(cubic_coeffs)
        ### We need to find the roots of 
        
        # Negative real values have angle pi, complex values have other, non-zero angle:
        real_pos_roots_mask = (np.angle(np_roots) == 0)

        real_roots = np_roots[real_pos_roots_mask].real

        minimum_roots = []
        for xi_root in real_roots:
            # Check that is minimum - second derivative should be positive at minima
            d2E_d2xi_at_root = self.lambified_d2e_dxi2(xi_root, Y_value, P_r_til_value)
            if d2E_d2xi_at_root > 0:
                minimum_roots.append(xi_root)
        
        
        
        if len(minimum_roots) != 1:
            raise ValueError("Failed to find real root - roots {}\n complex angles - {}".format(np_roots, np.angle(np_roots)))
        
        return minimum_roots[0]

    def generate_lambdified_grad_and_hessian_function(self):
        gradient_element_funcs = []
        for i, i_var in enumerate([self.alpha_til, self.beta_til]):
            grad_element_sym = self.E_til.diff(i_var)

            gradient_element_funcs.append(
                sympy.lambdify((self.alpha_til, self.beta_til,self.Y,self.P_r_til),grad_element_sym),
            )
        
        hessian_mat_element_funcs = []
        for i, i_var in enumerate([self.alpha_til, self.beta_til]):
            hessian_elfunc_row = []
            for j, j_var in enumerate([self.alpha_til, self.beta_til]):
                hessian_el_symbol = self.E_til.diff(i_var, j_var)
                # return hessian_el_symbol
                hessian_elfunc_row.append(
                    sympy.lambdify((self.alpha_til, self.beta_til,self.Y,self.P_r_til), hessian_el_symbol)
                )
            hessian_mat_element_funcs.append(hessian_elfunc_row)

        def gradient_function(alpha_t_value, beta_t_value, Y_value, P_r_til_value):
            grad_numeric = np.zeros((2,),dtype=float)
            for i in range(2):
                grad_numeric[i] = gradient_element_funcs[i](
                    alpha_t_value,
                    beta_t_value,
                    Y_value,
                    P_r_til_value,
                )
            return grad_numeric
            
        def hessian_function(alpha_t_value, beta_t_value, Y_value, P_r_til_value):
            hessian_numeric = np.zeros((2,2), dtype=float)
            for i in range(2):
                for j in range(2):
                    hessian_numeric[i][j] = hessian_mat_element_funcs[i][j](
                        alpha_t_value,
                        beta_t_value,
                        Y_value,
                        P_r_til_value,
                    )

            return hessian_numeric
        
        
        return gradient_function, hessian_function
    
    def find_rest_size_of_hexagon(
        self,
        A_0,
        P_0,
        K,
        Gamma,
    ):
        if A_0 <= 0:
            raise ValueError("Invalid value A_0 = {} - must be positive".format(A_0))
        if K <= 0:
            raise ValueError("Invalid value K = {} - must be positive".format(K))
        
        mu_val = np.sqrt(A_0) # Scaling factor to get tilded variables
        
        Y_val = Gamma / (K*A_0)
        P_r_til_val = P_0 / mu_val
        
        
        eq_vals = self.get_important_vals_at_givenpoint(Y_val, P_r_til_val)
        
        xi_til_eq = eq_vals["xi_til_eq"] # Equilibrium side length, in tilded variable
        poisson_ratio = eq_vals["poisson_ratio"] # Same in either variable
        
        side_length_eq = xi_til_eq * mu_val
        rest_area = (3*np.sqrt(3)/2) * (side_length_eq ** 2)
        
        # eigenvalues, eigenvectors = np.linalg.eig(eq_vals[""]
        if eq_vals["poisson_ratio"] > 1:
            raise Exception("Poisson ratio more than one, which means floppy regime - both A0 and P0 can be satisfied in floppy state...")
        
        return {
            "rest_side_length": side_length_eq,
            "rest_area": rest_area,
            "poisson_ratio": poisson_ratio,
        }

    def get_important_vals_at_givenpoint(
        self,
        Y_value, # Defined as Gamma/(K*A0).
        P_r_tilvalue, # Defined in sijie's paper as p_0=P0/sqrt(A0).
    ):
        xi_til_eq = self.find_xi_til_equilibrium_numeric(
            Y_value,
            P_r_tilvalue,
        )
        
        gradient_at_eq = self.get_gradient(
            xi_til_eq,
            xi_til_eq,
            Y_value,
            P_r_tilvalue,
        )
        
        hessian_at_eq = self.get_hessian(
            xi_til_eq,
            xi_til_eq,
            Y_value,
            P_r_tilvalue,
        )
        
        poisson_ratio = hessian_at_eq[0][1] / hessian_at_eq[1][1]
        
        return {
            "poisson_ratio": poisson_ratio,
            "xi_til_eq": xi_til_eq,
        }


    def plot_values_of_Y_and_Pr_til(self):
        Y_min = 0.001
        Y_max =  100
        Y_range = np.logspace(
            np.log10(Y_min),
            np.log10(Y_max),
            200,
        )
        
        Pr_t_min = -10
        Pr_t_max = 10
        Pr_t_range = np.linspace(
            Pr_t_min,
            Pr_t_max,
            200,
        )
        # Pr_t_range = np.logspace(
        #     np.log10(Pr_t_min),
        #     np.log10(Pr_t_max),
        #     100,
        # )
        Y_meshvals, Pr_t_meshvals = np.meshgrid(
            Y_range,
            Pr_t_range,
        )

        poisson_val_meshvals = np.zeros_like(Y_meshvals,dtype=float)
        for i in range(Y_meshvals.shape[0]):
            for j in range(Y_meshvals.shape[1]):
                try:
                    vals_at_pt = self.get_important_vals_at_givenpoint(
                        Y_meshvals[i][j],
                        Pr_t_meshvals[i][j],
                    )
                    poisson_val_meshvals[i][j] = vals_at_pt["poisson_ratio"]
                except Exception as e:
                    pass
                    poisson_val_meshvals[i][j]=np.nan
                    # print("Y: {}".format(Y_meshvals[i][j]))
                    # print("Pr_t: {}".format(Pr_t_meshvals[i][j]))
                    # pass
                    # print(str(e))
                    # break
                # except:
                #     pass


        ## Plot stuff
        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize=(4,4))

        p_ratio_max_magnitude = np.abs(poisson_val_meshvals[~np.isnan(poisson_val_meshvals)]).max()

        cmap = colormaps['bwr']
        cmap.set_bad(color="black")

        im_poisson_ratio = ax.pcolormesh(
            Y_meshvals,
            Pr_t_meshvals,
            poisson_val_meshvals,
            cmap=cmap ,
            vmin= -p_ratio_max_magnitude,
            vmax=  p_ratio_max_magnitude,
        )
        # mpl.colormaps['viridis']

        cbar_poisson  =fig.colorbar(im_poisson_ratio, orientation='vertical', label=r"Color: Poisson Ratio=$\nu$")#, format=LogFormatterMathtext())
        
        ax.set_xscale('log')
        # ax.set_yscale('log')
        
        ax.set_xlabel(r"$\Gamma/(KA_r)$")
        ax.set_ylabel(r"$p_0=P_r/\sqrt{A_r}$")

        p0_crit = np.sqrt(24/np.sqrt(3))
        
        for Y_val, col in zip([
            1/3.4641,
            1/34.641,
            1/346.31,
        ], [
            "salmon",
            "orange",
            "yellow",
        ]):
            ax.axvline(
                Y_val,
                c=col,
                linestyle='--',
                label=r"$\Gamma/(KA_r)={:2f}$".format(Y_val),
            )
            ax.scatter(x=Y_val,y=p0_crit, color='black')


        ax.axhline(p0_crit,c='black',linestyle='--',label=r"$p_0=p_c \approx 3.722$")
        
        plt.subplots_adjust(right=0.7)
        ax.legend(loc='upper left', bbox_to_anchor=(1.3, 1.0))
        
        plt.show()


    
    

# @classmethod
def find_regular_hexagonal_rests_area():
    sympy.init_printing()
    
    
    hex_m = RegHexagonalModel()
    print(sympy.Eq(
        sympy.Function('E')(hex_m.alpha, hex_m.beta, hex_m.A_r, hex_m.P_r, hex_m.K, hex_m.Gam),
        hex_m.energy_basic,
    ))

    print(hex_m.tilde_mu_substitutions)
    
    

    print(sympy.Eq(
        sympy.Function("E")(hex_m.alpha_til, hex_m.beta_til, hex_m.K_til, hex_m.Gam_til ,hex_m.P_r_til),
        hex_m.E_mu_subbed,
    ))
    

    print(sympy.Eq(
        sympy.Function(r"\tilde{E}")(hex_m.alpha_til, hex_m.beta_til, hex_m.Y,hex_m.P_r_til),
        hex_m.E_til,
    ))
    

    print(sympy.Eq(
        sympy.Function(r"\tilde{E}")(hex_m.xi_til, hex_m.Y,hex_m.P_r_til),
        hex_m.E_til_iso,
    ))


    print(sympy.Eq(
        sympy.Derivative(sympy.Function(r"\tilde{E}")(hex_m.xi_til), hex_m.xi_til),
        hex_m.dE_dxi_symexp
    ))
    
    print(hex_m.dE_dxi_polynomial)
    
    print(hex_m.d2E_dxi2_sym)
    
    
    

    print(hex_m.dE_dxi_coeff_functions)
    
    hex_m.plot_values_of_Y_and_Pr_til()
    
# %matplotlib widget


