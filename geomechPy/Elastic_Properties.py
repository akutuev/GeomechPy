import numpy as np
import scipy as sp

class Elastic_Properties:
	#Velocity2SLowness
	def slowness(vs, vp,unit_in,unit_o):
		if unit_in == 'm/s' :
			aa = 304800/vs
			bb = 304800/vp
		elif unit_in == 'km/s' and unit_o == 'us/m':
			
			aa = 304800/(vs*1000)
			bb = 304800/(vp*1000)
		return aa , bb


	#Slowness2Velocity
	def velocity(ds, dt,unit):
		if unit_in == 'us/ft':
			aa = 304800/ds
			bb = 304800/dt
		elif unit_in == 'km/s':
			aa = 304.800/ds
			bb = 304.800/dt
		return aa , bb
		
	def DynProp_Sonic(dt, ds, rho):
		G = rho*(304800/ds)**2/1000000
		M = rho*(304800/dt)**2/1000000
		vpvs = ds/dt
		
		
		K = M - 4/3*G
		E = G*(3*M  - 4*G)/(M - G)
		L = M- 2*G
		PR = (M - 2*G)/(2*M - 2*G)	
		
		return G, M, vpvs, E, K, PR , L

	def DynProp_Vel(vp, vs, rho,unit_out='GPA',digits='3'):
		round = int(digits)
		if unit_out == 'PSI':
			unit_multi =0.145038
		else:
			unit_multi = 1
		G = np.round((rho*(vs)**2/1e6)*unit_multi,round)
		M = np.round((rho*(vp)**2/1e6)*unit_multi,round)
		
		vpvs = np.round(vp/vs,round)
		K = np.round(M - 4/3*G,round)
		E = np.round(G*(3*M  - 4*G)/(M - G),round)
		L = np.round(M- 2*G,round)
		PR = np.round((M - 2*G)/(2*M - 2*G),round)

		return G, M, vpvs, E, K, PR , L

	def Elastic_Property_Conversion(aa, bb,input='E-PR',unit_out = 'GPa',digits='3'):
		"""
		Conversion between dynamic elastic properties.
		
		Args:
			a,b: Pair of Elastic property Young's Modulus - Y, Bulk Modulus - K, Poisson's Ratio - PR, Shear Modulus - G, Compressional Modulus - M , Lame Parameter - L
			Input can have any units as long as consistent: Y,K,G,M,L --> Pressure Units (psi, GPa, bar etc) PR--> Unitless
			input: Define what the combination of input elasti properties is: Default 'E-PR', Other options: 'K-E','K-L','K-G','K-PR','K-M'
			digits: Round the output to number of digits Default: '3'
		
		
		Returns:
			Output are the four elastic properties different from the input pair
			Example: if inpt is 'E-PR' output will be 'K-G-L-M'
		Reference:
		https://en.wikipedia.org/wiki/Elastic_modulus
		"""

		a = np.array(aa)
		b = np.array(bb)
		if unit_out == ():
			unit_multi = 1
		elif unit_out == 'Mpsi':
			unit_multi =145038/1e6
		elif unit_out == 'psi':
			unit_multi =145038
		elif unit_out == 'bar':
			unit_multi = 10000						
		elif unit_out == 'Pa':
			unit_multi = 1000000

		round = int(digits)
		if input == 'K-E':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			L  = 3*a*(3*a - b)/(9*a - b)
			G  = 3*a*b/(9*a-b)
			PR = (3*a - b)/(6*a)
			M  = 3*a*(3*a + b)/(9*a - b)
			return L, G, PR, M
		elif input == 'K-L':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			E  = 9*a*(3*a-b)/(3*a-b)
			G  = 3*(a-b)/2
			PR = b/(3*a-b)
			M  = 3*a - 2*b
			return E, G, PR, M
		elif input == 'K-G':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			E  = 9*a*b/(3*a+b)
			L  = a - 2*b/3
			PR = (3*a - 2*b)/(2*(3*a+b))
			M  = a +4*b/3
			return E, L, PR, M
		elif input == 'K-PR':
			a = np.array(aa)*unit_multi
			b = np.array(bb)
			E  = 3*a*(1-2*b)
			L  = 3*a*b/(1+b)
			G = 3*a*(1-2*b)/(2*(1+b))
			M  = 3*a*(1-b)/(1+b)
			return E, L, G, M
		elif input == 'K-M':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			E  = 9*a*(b-a)/(3*a+b)
			L  = (3*a-b)/2
			G = (3*(b-a))/4
			PR  = (3*a-b)/(3*a+b)
			return E, L, G, PR
		elif input == 'E-L':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			R = sqrt(a**2 +9*b**2 +2*ab)
			K  = (a+3*b+R)/6
			G  = (a-3*b+R)/4
			PR = 2*b/(a+b+R)
			M  = (a-b+R)/2
			return K, G, PR, M
		elif input == 'E-G':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			K  = a*b/(3*(3*b-a))
			L  = (b*(a-2*b))/(3*b-a)
			PR = 0.5*(a/b) - 1
			M  =  (b*(4*b-a))/(3*b-a)
			return K, L,PR , M 
		elif input == 'E-PR':
			a = np.array(aa)*unit_multi
			b = np.array(bb)
			K = a/(3*(1-2*b))
			L = a*b/((1+b)*(1-2*b))
			G = a/(2*(1+b))
			M =  (a*(1-b))/((1+b)*(1-2*b))
			return K, L, G, M
		elif input == 'E-M':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			S = np.sqrt(a**2 + 9*b**2 -10*a*b)
			K  = (3*b - a + S)/6
			L  = (b-a+S)/4
			G = (3*b+a-S)/8
			PR = (a-b+S)/(4*b)
			return K,L ,G ,PR 
		elif input == 'L-G':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			K  = a + (2/3)*b
			E  = b*(3*L +2*b)/(a+b)
			PR = 0.5*a/(a-b)
			M  = a +2*b
			return K,E , PR, M 
		elif input == 'L-PR':
			a = np.array(aa)*unit_multi
			b = np.array(bb)
			K  = a*(1+b)/(3*b)
			E  = (a*(1+b)*(1-2*b))/b
			G  = (b*(1-2*b))/(2*b)
			M  = (a*(1-b))/b 
			return K, E,G , PR
		elif input == 'L-M':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			K  = (b+2*a)/3
			E  = ((b-a)*(b+2*a))/(b+a)
			G =  (b-a)/2
			PR = a/(b+a) 
			return K,E ,G ,L 
		elif input == 'G-PR':
			a = np.array(aa)*unit_multi
			b = np.array(bb)
			K  = (2*a*(1+b))/(3*(1-2*b))
			E  = 2*a*(1+b)
			L = 2*a*b/(1-2*b)
			M = (2*a*(1-b))/(1-2*b)
			return K,E ,L ,M 
		elif input == 'G-M':
			a = np.array(aa)*unit_multi
			b = np.array(bb)*unit_multi
			K  = b- (4/3)*a
			E  = (a*(3*b-4*a))/(b-a)
			L  = b -2*a 
			PR  = (b -2*a)/(2*b -2*a)
			return K,E,L,M 
		elif input == 'PR-M':
			a = np.array(aa)
			b = np.array(bb)*unit_multi
			K  = (b*(1+a))/(3*(1-a))
			E  = (b*(1+a)*(1-2*a))/(1-a)
			L = b*a/(1-a)
			G = (b*(1-2*a))/(2*(1-a))
			return K,E,L,G 
		

	def UCS_correlations(aa,Method='Horsrud',Input='DTCO',unit_in='usft',unit_out='MPa',digits='1'):
		round = int(digits)
		if unit_out == 'PSI':
			unit_multi =145.038
		else:
			unit_multi = 1
		if Input == 'DTCO':
			if Method == 'Horsrud':
				UCS = np.round((0.77*pow(304.8/aa,2.93))*unit_multi,round) #Horsrud
			elif Method == 'McNally':
				UCS = np.round((1200.0*np.exp(-0.036*aa))*unit_multi,round) #McNally
			else:
				print("something must be wrong")
		elif Input == 'YME_STA':
			if Method == 'Plumb_YME_STA':
				UCS = np.round((4.424*aa)*unit_multi,round) #McNally
		else:
			print("something must be wrong")
			return UCS
		
		


