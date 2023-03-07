    # Fraction of energy between two frequencies IN OBSERVER FRAME
    def energy_frac(self, lowerlim, upperlim):
        p = 2.2
        nu_c, nu_m = self.nu_c, self.nu_m

        # Convert frequencies to comoving frame
        doppler_factor = np.exp(np.arccosh(self.gamma_r))
        lowerlim /= doppler_factor
        upperlim /= doppler_factor
        
        # Normalization constants in regions B and C
        normB = nu_c**(5/6)
        normBslow = nu_m**(p/2 - 1/6)
        normC = normB * nu_m**((p-1)/2)
        normCslow = normBslow*nu_c**0.5

        # Domain is divided in regions A, B, C
        # There are six possibilities
        def aa(lowerlim, upperlim):
            return 3/4 * (upperlim**(4/3) - lowerlim**(4/3))
        
        def bb(lowerlim, upperlim):
            if nu_c < nu_m: # fast cooling
                return 2*normB * (upperlim**(1/2) - lowerlim**(1/2))
            else:           # slow cooling
                return 2*normBslow/(3-p) * (upperlim**((3-p)/2) - lowerlim**((3-p)/2))
        
        def cc(lowerlim, upperlim):
            if nu_c < nu_m:
                return 2*normC/(2-p) * (upperlim**(1-p/2) - lowerlim**(1-p/2))
            else:
                return 2*normCslow/(p-2) * (lowerlim**(-p/2+1) - upperlim**(-p/2+1))
        
        def ab(lowerlim, upperlim):
            if nu_c < nu_m:
                return aa(lowerlim, nu_c) + bb(nu_c, upperlim)
            else:
                return aa(lowerlim, nu_m) + bb(nu_m, upperlim)
        
        def ac(lowerlim, upperlim):
            if nu_c < nu_m:
                return aa(lowerlim, nu_c) + bb(nu_c, nu_m) + cc(nu_m, upperlim)
            else:
                return aa(lowerlim, nu_m) + bb(nu_m, nu_c) + cc(nu_c, upperlim)
        
        def bc(lowerlim, upperlim):
            if nu_c < nu_m:
                return bb(lowerlim, nu_m) + cc(nu_m, upperlim)
            else:
                return bb(lowerlim, nu_c) + cc(nu_c, upperlim)
        
        if nu_c < nu_m:
            if lowerlim < nu_c:
                if upperlim < nu_c:
                    f = aa
                elif upperlim < nu_m:
                    f = ab
                else:
                    f = ac
            elif lowerlim < nu_m:
                if upperlim < nu_m:
                    f = bb
                else:
                    f = bc
            else:
                f = cc
        else:
            if lowerlim < nu_m:
                if upperlim < nu_m:
                    f = aa
                elif upperlim < nu_c:
                    f = ab
                else:
                    f = ac
            elif lowerlim < nu_c:
                if upperlim < nu_c:
                    f = bb
                else:
                    f = bc
            else:
                f = cc

        return f(lowerlim, upperlim) / ac(0, np.infty)

    
    def energy_between(self, lowerlim, upperlim):
        return self.energy * self.energy_frac(lowerlim, upperlim)