# Quantum algorithm for linear approximation of vectorial Boolean function

The C and python code implements the Proposed Quantum algorithm for linear approximation of a non-linear vectorial Boolean function.
<a href=https://quest.qtechtheory.org/>QuEST</a> library is used to simulate quantum environment.<br/>
<a href=https://www.ibm.com/quantum-computing/>IBM's Quantum Computer</a> is used for executing the proposed quantum circuit.<br/><br/>

"q_linapx.c" contains the C code of the proposed algorithm using QuEST library. "q_linapx.exe" is the respective executable for windows.<br/>
To run "q_linapx.exe", simply start the command prompt and traverse to the directory that contains "q_linapx.exe" and type "q_linapx.exe" (without quotes) and press enter to get a high probability linear approximation. The code prints a high probability linear approximation (in hexadecimal and decimal form).<br/>

"q_linapx.c" contains the C code of the proposed algorithm using QuEST library which repeats the proposed quantum circuit for a number of times and measure the result. "q_linapx.exe" is the corresponding executable for windows.<br/>
Usage of "q_linapx.exe": Type in command prompt "q_linapx.exe <times>" without quotes. Replace <times> by a positive integer to repeat the quantum circuit <times> times.<br/>

"q_linapx.ipynb" contains the python code of the proposed circuit to run on IBM's real quantum computer. Use jupyter notebook to run the code. Note that token ID from IBM is required in "save_account()" function inside the file.<br/>
