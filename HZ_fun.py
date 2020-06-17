def HZ_in_out_2013(T_eff,Lum):
	"""
	
	HABITABLE ZONES AROUND MAIN-SEQUENCE STARS ! 
	
	Parameters:
	-----------

	T_eff: float;
	    Effective temperatures
	    Unit [Sun]

	Lum: float;
	    Luminosity
	    Unit: [Sun]

	Returns:
	--------
	https://arxiv.org/pdf/1301.6674.pdf: 

	Recent Venus, Runaway Greenhouse, Moist Greenhouse, Maximum Greenhouse,\
	 Early Mars , conservative[Moist Greenhouse, Maximum Greenhouse] ,\
	  optimistic[Recent Venus, Early Mars]

	"""
	T_star = T_eff - 5780
	S1, S2, S3, S4, S5 = 1.7753, 1.0512, 1.0140, 0.3438, 0.3179
	a1, a2, a3, a4, a5  = 1.4316E-4, 1.3242E-4, 8.1774E-5, 5.8942E-5, 5.4513E-5
	b1, b2, b3, b4, b5  = 2.9875E-9, 1.5418E-8, 1.7063E-9, 1.6558E-9, 1.5313E-9
	c1, c2, c3, c4, c5  = -7.5702E-12, -7.9895E-12, -4.3241E-12, -3.0045E-12, -2.7786E-12
	d1, d2, d3, d4, d5  = -1.1635E-15, -1.8328E-15, -6.6462E-16, -5.2983E-16, -4.8997E-16
	s_eff1 = S1 + a1*T_star + b1*T_star**2 + c1*T_star**3 + d1*T_star**4
	s_eff2 = S2 + a2*T_star + b2*T_star**2 + c2*T_star**3 + d2*T_star**4
	s_eff3 = S3 + a3*T_star + b3*T_star**2 + c3*T_star**3 + d3*T_star**4
	s_eff4 = S4 + a4*T_star + b4*T_star**2 + c4*T_star**3 + d4*T_star**4
	s_eff5 = S5 + a5*T_star + b5*T_star**2 + c5*T_star**3 + d5*T_star**4
	d1, d2 ,d3, d4, d5 = (Lum/s_eff1)**(0.5), (Lum/s_eff2)**(0.5), (Lum/s_eff3)**(0.5),\
	 (Lum/s_eff4)**(0.5), (Lum/s_eff5)**(0.5)
	conservative = [d3, d4]
	optimistic = [d1, d5]
	return d1, d2, d3, d4, d5, conservative, optimistic
