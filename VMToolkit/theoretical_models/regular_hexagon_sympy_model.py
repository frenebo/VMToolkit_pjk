import sympy
import numpy as np


"""
Everything here is adapted from the code in this notebook:
2024_11_27_sijie_system.ipynb
which repeated the previous caculations I did in terms of the energy paraameteers used by Sijie in his paper,
and by most people working with vertex models.

This model finds the equilibrium size of a regular hexagon cell, then finds shear/bulk moduli, and other elastic
properties in the linear regime (e.g. small deformations in x and y.)

Note: this model ONLY allows scaling along x and y axis - so non-affine deformations are NOT allowed.
A more complete model would allow vertices more freedom to deform, beyond infinitesimal affine deformations
from the regular hexagon it starts at.
"""
class TheoreticalRegularHexModel:
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
    
    # @NOTE - supposed to be evaluated at equilibrium xi value!
    @classmethod
    def create_bulk_modulus_basic(
        cls,
        energy_regular_exp,
        alpha,
        beta,
        xi_equilibrium,
        A_r,
        P_r,
        K,
        Gam
    ):
        xi = xi_equilibrium
        # Energy, if alpha and beta scaling parameters are set equal to xi - isotropic shrink or stretch
        energy_regular_iso = energy_regular_exp.subs([
            (alpha, xi),
            (beta, xi),
        ])
        
        
        # pressure = dE/dA 
        # dA = Perimeter * d(center-to-edge distance) = 6*xi*d(xi*sqrt(3)/2) = 3*sqrt(3)*xi*dxi
        # pressure = dE/(3*sqrt(3)*xi*dxi) = 1/(3*sqrt(3)*xi)*dE/dxi
        
        # We need dE/dgamma to find pressure
        dE_regular_iso_dxi_symexp = energy_regular_iso.diff(xi).expand()
        
        pressure_reg_iso_exp = ( (-1) / (3 * sympy.sqrt(3) * xi) ) * dE_regular_iso_dxi_symexp
        pressure_reg_iso_exp = pressure_reg_iso_exp.simplify()
        
        # Bulk modulus is defined as:
        # B = -A d(pressure)/dA
        # B = -A d(pressure)/dxi / (dA/dxi)
        # Fining dA/dxi:
        # A = (3*sqrt(3)/2)*xi^2
        A_exp = ((3*xi*xi)*sympy.sqrt(3))/2
        dA_dxi_exp = A_exp.diff(xi)
        dpressure_dxi_reg_iso_exp = pressure_reg_iso_exp.diff(xi).simplify()
        
        B_bulkmod_reg_iso = (-1) * A_exp * dpressure_dxi_reg_iso_exp / dA_dxi_exp
        
        return B_bulkmod_reg_iso
        
    
    def __init__(self):
        alpha, beta = sympy.symbols("alpha beta", real=True, positive=True)
        A_r, P_r = sympy.symbols("A_r P_r", real=True, positive=True)
        K, Gam = sympy.symbols("K Gamma", real=True)
        mu = sympy.symbols("mu", real=True, positive=True)
        
        Y = sympy.symbols("Y", real=True, positive=True)
        
        self.alpha = alpha
        self.beta = beta
        self.A_r = A_r
        self.P_r = P_r
        self.K = K
        self.Gam = Gam
        self.mu = mu
        self.Y = Y
        # self.xixi_equililbrium_sym = xi_equililbrium_sym
        
        self.energy_basic = self.create_energy_symexp(
            alpha,beta, A_r, P_r, K, Gam,
        )
        
        # Tilde variables make the math easier

        Gam_til, K_til, alpha_til =  sympy.symbols(
                                                        r"\tilde{\Gamma} \tilde{K} \tilde{\alpha}",
                                                        real=True,
                                                        positive=True)
                                                        
        beta_til, P_r_til, A_r_til =  sympy.symbols(
                                                        r"\tilde{\beta}  \tilde{P_r} \tilde{A_r}",
                                                        real=True,
                                                        positive=True)
        xi_til =  sympy.symbols(r"\tilde{\xi}", real=True, positive=True)

        self.Gam_til = Gam_til
        self.K_til = K_til
        self.alpha_til = alpha_til
        self.beta_til = beta_til
        self.P_r_til = P_r_til
        self.A_r_til = A_r_til
        self.xi_til = xi_til
        
                                                        
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
        
        self.get_bulk_modulus_numerical = self.generate_bulk_mod_numeric_function()
        
        
    def generate_bulk_mod_numeric_function(self):
        xi_equilibrium_sym =  sympy.symbols(r"x_{eq}", real=True, positive=True)

        
        # @NOTE - supposed to be evaluated at equilibrium xi value!
        bulk_modulus_basic_exp = self.create_bulk_modulus_basic(
            energy_regular_exp=self.energy_basic,
            alpha=self.alpha,
            beta=self.beta,
            xi_equilibrium=xi_equilibrium_sym,
            A_r=self.A_r,
            P_r=self.P_r,
            K=self.K,
            Gam=self.Gam
        )
        print("Bulk modulus:")
        print(bulk_modulus_basic_exp)
        lambdified_reg_vars = sympy.lambdify(
            (
                xi_equilibrium_sym,
                self.A_r,
                self.P_r,
                self.K,
                self.Gam,
            ),
            bulk_modulus_basic_exp
        )
        
        def numeric_func(
            xi_reg_eq_value,
            
            gamma_reg_value,
            K_reg_value,
            
            P_r_reg_value, # Defined in sijie's paper as p_0=P0/sqrt(A0).
            A_r_reg_value,
        ):
            # print("Xi: {}".format(xi_reg_eq_value))
            # print("gamma: {}".format(gamma_reg_value))
            # print("K: {}".format(K_reg_value))
            
            # print("P_r: {}".format(P_r_reg_value))
            # print("A_r: {}".format(A_r_reg_value))
            return lambdified_reg_vars(
                xi_reg_eq_value,
                A_r_reg_value,
                P_r_reg_value,
                K_reg_value,
                gamma_reg_value,
            )
            
        return numeric_func
        
        
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
    
    def find_elastic_props_of_hexagon(
        self,
        A_0_num,
        P_0_num,
        K_num,
        gamma_num,
    ):
        if A_0_num <= 0:
            raise ValueError("Invalid value A_0_num = {} - must be positive".format(A_0_num))
        if K_num <= 0:
            raise ValueError("Invalid value K_num = {} - must be positive".format(K_num))
        
        mu_val = np.sqrt(A_0_num) # Scaling factor to get tilded variables
        
        
        eq_vals = self.get_important_vals_at_givenpoint(
            gamma_reg_value=gamma_num,
            P_r_reg_value=P_0_num,
            K_reg_value=K_num,
            A_r_reg_value=A_0_num,
        )
        
        xi_til_eq     = eq_vals["xi_til_eq"] # Equilibrium side length, in tilded variable
        poisson_ratio = eq_vals["poisson_ratio"] # Same in either variable
        bulk_modulus = eq_vals["bulk_modulus"]
        
        side_length_eq = xi_til_eq * mu_val
        rest_area = (3*np.sqrt(3)/2) * (side_length_eq ** 2)
        
        if poisson_ratio <= -1:
            raise ValueError("Poisson ratio less than or equal to minus one!")
            
        if eq_vals["poisson_ratio"] > 1:
            print("POISSON RATIO GREATER THAN ONE!")
            rest_area = A_0_num
            side_length_eq = np.sqrt(rest_area  * 2/(3*np.sqrt(3)) )
            
            poisson_ratio = 1.0
            
            shear_modulus = 0.0
        else:
            # This is how it works in 3d:
            # G_3d = shear modulus, B=bulk modulus
            # 2G(1+nu)=3B(1-2nu)
            #      3B(1-2nu)
            # G_3d = -----------
            #       (2+2nu)
            # But in 2d...
            # G = B*(1-nu)/(1+nu)
            # Can find by redoing derivation in landau lifshitz elasticity Ch1.Sec4, but in 2D.
            # So (1/3)delta_ij => 1/2 delta_ij, take sums just over x&y, etc.
            shear_modulus = bulk_modulus*(1 - poisson_ratio) / (1 + poisson_ratio)
        # if sshear
        
        print("Bulk modulus: {}".format(bulk_modulus))
        print("Shear modulus: {}".format(shear_modulus))
        print("Poisson ratio: {}".format(poisson_ratio))
        
        return {
            "rest_side_length": side_length_eq,
            "rest_area": rest_area,
            "poisson_ratio": poisson_ratio,
            "bulk_modulus": bulk_modulus,
            "shear_modulus": shear_modulus,
        }

    def get_important_vals_at_givenpoint(
        self,
        gamma_reg_value,
        P_r_reg_value,
        K_reg_value,
        A_r_reg_value,
    ):
        mu_val = np.sqrt(A_r_reg_value)
        # Finding equilibrium and poisson ratio is easier in tilded coordinates
        Y_value = gamma_reg_value / (K_reg_value*A_r_reg_value)
        P_r_tilvalue = P_r_reg_value / mu_val
        
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
        
        # Change of variables back to normal, to calculate bulk modulus.
        # Bulk moduus is not dimensionless, so we need to change back to do calculations
        xi_nontilded_eq = xi_til_eq * mu_val
        
        bulk_modulus = self.get_bulk_modulus_numerical(
            xi_reg_eq_value=xi_nontilded_eq,
            
            A_r_reg_value=A_r_reg_value,
            P_r_reg_value=P_r_reg_value,
            
            gamma_reg_value=gamma_reg_value,
            K_reg_value=K_reg_value,
        )
        
        
        poisson_ratio = hessian_at_eq[0][1] / hessian_at_eq[1][1]
        
        return {
            "poisson_ratio": poisson_ratio,
            "xi_til_eq": xi_til_eq,
            "bulk_modulus": bulk_modulus
        }


    def plot_values_of_Y_and_Pr_til(self):
        from matplotlib import colormaps
        from matplotlib import pyplot as plt
        
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


    
    

def find_regular_hexagonal_rests_area():
    print("This function hasn't been tested in a while...")
    sympy.init_printing()
    
    
    hex_m = TheoreticalRegularHexModel()
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
    