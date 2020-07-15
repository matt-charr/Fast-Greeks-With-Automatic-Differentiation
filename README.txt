************************************************************************************************
*COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI)*
*ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM******
*OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM********
************************************************************************************************

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-> PROJECT NAME


	Fast Greeks Computation with Automatic Differentiation Tools


-> REQUIREMENT AND INSTALLATION


	- C++17 compiler
	- Usage of 'llapack', 'blas' and 'wall'
	- 'cmake' tool version 3.0.0 at least


-> EXECUTABLES FILES


	Automatic Differentiation exhibition:

		- 'demoTangentMode': Tangent Mode of AD demo
		- 'demoAdjointMode': Adjoint Mode of AD demo

	Processing time save:

		- 'basketOptionTime': Basket Option greeks calculation
		- 'bestOfAsianOptionTime': Best Of Asian Option greeks calculation
		- 'correlationTime': Correlation greeks of Basket Option calculation	
		
	Results consistency test:

		- 'bestOfAsianOptionShow': Monte Carlo convergence of Best Of Asian Option greeks calculation
		- 'basketOptionShow': Monte Carlo convergence of Basket Option greeks calculation
		- 'correlationShow': Monte Carlo convergence of Basket Option correlation greeks calculation


 -> WHERE TO START FOR A DEMO ?


	- The user can first test the efficiency of tangent mode with 'demoTangentMode.cpp':
		
		. Start by writing the function to differentiate (scalar valued functions with multiple variables in input are required)
		. The user can use any binary/unary operators and external functions he want -- PROVIDED the operation/function was overloaded in the tangent mode library
		  [To overload an operation or a function, refer to 'tgt_double.hpp']
		. The user is asked the value of the inputs and the coefficients associated to each input variable in the linear combination of derivatives (any scalars)
		. In order to introduce a new input variable, do not forget to add the corresponding code lines.
		. To declare a new 'tgt_double' object, use the command 'X = make0<tgt_double> (X_value, X_coefficient)'

	- The user can then test the efficiency of adjoint mode with 'demoAdjointMode.cpp':
		
		. Start by writing the function to differentiate (scalar valued functions with multiple variables in input are required)
		. The user can use any binary/unary operators and external functions he want -- PROVIDED the operation/function was overloaded in the adjoint mode.
		[To overload an operation or a function, refer to 'adj_double.hpp']
		. The user is asked the value of the inputs and the seed (proportionality coefficient in the derivative)
		. In order to introduce a new input variable, do not forget to add the corresponding code lines.
		. To declare a new 'adj_double' object, use the command 'X = make0<adj_double> (X_value, X_coefficient)'

	- The user can assess the processing time to compute basket option/best of asian option prices and greeks using the three methods exposed (Finite Difference Method, 
	  Tangent Mode of AD, Adjoint Mode of AD) with the files 'basketOptionTime.cpp', 'correlationTime.cpp' and 'bestOfAsianOptionTime.cpp'

	- The user can assess the coherence of the result with basket option/best of asian option prices and greeks using the three methods exposed (Finite Difference Method, 
	  Tangent Mode of AD, Adjoint Mode of AD) with the files 'basketOptionShow.cpp', 'correlationShow.cpp' and 'bestOfAsianOptionShow.cpp'


--> CONTRIBUTING


	- The user is free to overload any operations/functions for his personal need by adding the coding lines in the header files 'adj_double.hpp' and 'std_double.hpp'
	  while keeping the coherence with the structure.

	- The user can add different payoffs following the structure presenting in the header files 'basketOption.hpp' and 'bestOfAsianOption.hpp'


-> COMMAND TO BUILD THE PROJECT


	In the directory where the project lies, do in the command line:

	- cd build
	- cmake ..
	- make


--> COMMAND TO EXECUTE A FILE


	- cd ../bin/file_name


-> TECHNOGLOGIES:

	
	- Armadillo (http://arma.sourceforge.net/docs.html)
		
		. Library to operate on matrix, vector, cube and field object

	- Mersenne lister generator (http://www.cplusplus.com/reference/random/mt19937/)

		. Library to produce pseudo random numbers with Mersenne Twister Generator


--> CREDITS:


	. Monte Carlo, random variable and Black Scholes classes are high inspired by the one built during the training session of Vincent Lemaire 
	(https://www.lpsm.paris/pageperso/lemaire/docs/restreint/MC_CPP.html - private access)


-> AUTHOR:


	Matthieu Charrier
	Student at MSc. Probability & Finance
	Sorbonne Université
	matthieucharrier1994@gmail.com

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------