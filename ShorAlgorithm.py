#!/usr/bin/env python
# coding: utf-8

# In[1]:


class ShorN6X5:
        
    def __init__(self, n=6, x=5, n_count=4, shots=1000, backend_name="ibmq_qasm_simulator", IsTwirl=False):        
        from math import pi
        
        # load IBMQ Accounts
        backend_use = BackendInit(backend_name)
        
        # variables
        self.n = n # factoring number
        self.x = x # modular
        self.n_count = n_count  # number of counting qubits
        self.shots = shots # number of shots
        self.backend_name = backend_use # simulation backend
        self.IsTwirl = IsTwirl # measurement error mitigation (twriling)
        self.Unitary_Twirl = {
            1:[0, 0, 0],
            2:[pi, 3*pi/2, pi/2],
            3:[pi, 0, pi],
            4:[0, pi/2, pi/2],
            5:[pi/2, 0, pi/2],
            6:[pi/2, pi, pi/2],
            7:[pi/2, pi/2, pi],
            8:[pi/2, pi/2, 0],
            9:[pi/2, 3*pi/2, pi],
            10:[pi/2, pi/2, 0],
            11:[pi/2, 0, 3*pi/2],
            12:[pi/2, pi, 3*pi/2]} # twirling parameters

    def ConstructQC(self, unitary_index=1):
        from qiskit import QuantumCircuit
        
        # construct quantum circuit
        qc = QuantumCircuit(self.n_count+3, self.n_count)
        
        for q in range(self.n_count):
            qc.h(q) # Initialize counting qubits in state |+>
        qc.x(self.n_count+2) # And auxiliary register in state |1>
        
        for q in range(self.n_count): # Do controlled-U operations
            qc.append(CxMod6(self.x, 2**q), [q] + [i+self.n_count for i in range(3)])
        qc.append(QftDagger(self.n_count), range(self.n_count)) # Do inverse-QFT
        
        if self.IsTwirl == True: # twirling
            qc.append(Twirling(self.Unitary_Twirl, unitary_index, self.n_count), [i for i in range(self.n_count)])

        qc.measure(range(self.n_count), range(self.n_count)) # measuremeN6X5nt

        return qc
        
    def SimulationQC(self, unitary_index_set=[1,5,9], eta=0.1):
        from qiskit import transpile, assemble
        from qiskit.visualization import plot_histogram
        
        if self.IsTwirl == False: # no twirling
            print('==================================================')
            print('Name of Backend :', self.backend_name) # Setting Simulation Backend
            print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number
            
            qc = self.ConstructQC() # construct quantum circuit
            backend = self.backend_name # setting backend
            t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile

            print("Transpiled Quantum Circuit Depth :", t_qc.depth())
            results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result
            
            # Counts Result
            counts = results.get_counts() # get counts
            print('Simulation Result (counts) :', counts)
            print('==================================================')
        
        else: # twirling
            ## Quantum pre-processing
            # Variable setting
            counts = {}
            
            # Apply 3 unitray indexes to twriling
            for unitary_index in unitary_index_set:
                print('==================================================')
                print('Unitary Index :', unitary_index) # Setting Unitary Index
                print('Name of Backend :', self.backend_name) # Setting Simulation Backend
                print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number

                qc = self.ConstructQC(unitary_index) # construct quantum circuit
                backend = self.backend_name # setting backend
                t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile
                
                print('Transpiled Quantum Circuit Depth :', t_qc.depth())
                results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result

                # Counts Result
                temp_counts = results.get_counts()
                print('Simulation Result (counts) :', temp_counts)
                
                counts = SumDictionary(counts, temp_counts)
                print('==================================================\n')

            ## Classical post-processing
            # Variable setting
            measurement_counts = []
            for i in range(self.n_count):
                measurement_counts.append([0, 0])
            counts_pp = {}
            
            for key in counts: # average count results
                counts[key] = round(counts[key]/3)
                
            # Calculate counts of [x_n | n]
            for i in range(self.n_count):
                for key in counts:
                    if key[i] == "0":
                        measurement_counts[i][0] += counts[key]
                    else:
                        measurement_counts[i][1] += counts[key]
                                    
            # Recover counts
            for key in counts:
                counts_pp[key] = counts[key]
                for i in range(self.n_count):
                    counts_pp[key] = counts_pp[key] * ( (1 - (eta*self.shots)/(2*measurement_counts[i][int(key[i])])) / (1-eta) )
                counts_pp[key] = round(counts_pp[key])
            counts = counts_pp
        
        return plot_histogram(counts)
        
    def QpexMod6(self):
        from qiskit import transpile, assemble
        from qiskit.transpiler import Layout
        
        qc = self.ConstructQC() # construct quantum circuit
        backend = self.backend_name  # setting backend        
        t_qc = transpile(qc, backend, optimization_level=2)  # quantum circuit transpile
        result = backend.run(assemble(t_qc, shots=1, memory=1), memory=True).result() # get result
        
        readings = result.get_memory()
        print("Register Reading: " + readings[0])
        
        phase = int(readings[0],2)/(2**self.n_count)
        print("Corresponding Phase: %f" % phase)
        
        return phase

    def RunShorAlgorithm(self):
        from fractions import Fraction
        from math import gcd
            
        factor_found = False
        attempt = 0
        while not factor_found:
            attempt += 1
            print("======================================================")
            
            phase = self.QpexMod6() # Phase = s/r
            print("\nAttempt %i:" % attempt)
            
            frac = Fraction(phase).limit_denominator(self.n) # Denominator should (hopefully!) tell us r
            r = frac.denominator
            print("Result: r = %i" % r)
            
            if phase != 0:
                guesses = [gcd(self.x**(r//2)-1, self.n), gcd(self.x**(r//2)+1, self.n)]
                print("Guessed Factors: %i and %i" % (guesses[0], guesses[1]))
                for guess in guesses:
                    if guess not in [1,self.n] and (self.n % guess) == 0: # Check to see if guess is a factor
                        print("*** Non-trivial factor found: %i ***" % guess)
                        factor_found = True        


# In[2]:


class ShorN21X4:
    
    def __init__(self, n=21, x=4, control_q=3, shots=1000, ancilla_q=2, backend_name="ibmq_qasm_simulator", IsTwirl=False):        
        from math import pi
        
        # load IBMQ Accounts
        backend_use = BackendInit(backend_name)

        # variables
        self.n = n # factoring number
        self.x = x # modular
        self.control_q = control_q # number of control qubits
        self.ancilla_q = ancilla_q # number of ancilla qubits
        self.shots = shots # number of shots
        self.backend_name = backend_use # simulation backend
        self.IsTwirl = IsTwirl # measurement error mitigation (twriling)
        self.Unitary_Twirl = {
            1:[0, 0, 0],
            2:[pi, 3*pi/2, pi/2],
            3:[pi, 0, pi],
            4:[0, pi/2, pi/2],
            5:[pi/2, 0, pi/2],
            6:[pi/2, pi, pi/2],
            7:[pi/2, pi/2, pi],
            8:[pi/2, pi/2, 0],
            9:[pi/2, 3*pi/2, pi],
            10:[pi/2, pi/2, 0],
            11:[pi/2, 0, 3*pi/2],
            12:[pi/2, pi, 3*pi/2]} # twirling parameters
    
    def ConstructQC(self, unitary_index=1):
        from qiskit import QuantumCircuit
        
        qc = QuantumCircuit(self.control_q+self.ancilla_q, self.control_q) # Circuit initialization
        for q in range(self.control_q): # Initialize control qubits
            qc.h(q)

        # Do controlled-U operations
        # Controlled U^1
        qc.cx(2, 4)

        # Controlled U^2
        qc.cx(1, 4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (1, 3, 4))
        qc.cx(4, 3)

        # Controlled U^4
        qc.x(4)
        qc.append(Margolous(self.control_q), (0, 4, 3))
        qc.x(4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (0, 3, 4))
        qc.cx(4, 3)

        # Do inverse-QFT
        qc.append(QftDagger(self.control_q), range(self.control_q))

        if self.IsTwirl == True: # twirling
            qc.append(Twirling(self.Unitary_Twirl, unitary_index, self.control_q), [i for i in range(self.control_q)])
        
        # Measure circuit
        qc.measure(range(self.control_q), range(self.control_q))
        
        return qc
    
    def SimulationQC(self, unitary_index_set=[1,5,9], eta=0.1):
        from qiskit import transpile, assemble
        from qiskit.visualization import plot_histogram
        
        if self.IsTwirl == False: # no twirling
            print('==================================================')
            print('Name of Backend :', self.backend_name) # Setting Simulation Backend
            print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number
            
            qc = self.ConstructQC() # construct quantum circuit
            backend = self.backend_name # setting backend
            t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile

            print("Transpiled Quantum Circuit Depth :", t_qc.depth())
            results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result
            
            # Counts Result
            counts = results.get_counts() # get counts
            print('Simulation Result (counts) :', counts)
            print('==================================================')

        else: # twirling
            ## Quantum pre-processing
            # Variable setting
            counts = {}
            
            # Apply 3 unitray indexes to twriling
            for unitary_index in unitary_index_set:
                print('==================================================')
                print('Unitary Index :', unitary_index) # Setting Unitary Index
                print('Name of Backend :', self.backend_name) # Setting Simulation Backend
                print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number

                qc = self.ConstructQC(unitary_index) # construct quantum circuit
                backend = self.backend_name # setting backend
                t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile
                
                print('Transpiled Quantum Circuit Depth :', t_qc.depth())
                results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result

                # Counts Result
                temp_counts = results.get_counts()
                print('Simulation Result (counts) :', temp_counts)
                
                counts = SumDictionary(counts, temp_counts)
                print('==================================================\n')

            ## Classical post-processing
            # Variable setting
            measurement_counts = []
            for i in range(self.control_q):
                measurement_counts.append([0, 0])
            counts_pp = {}
            
            for key in counts: # average count results
                counts[key] = round(counts[key]/3)
                
            # Calculate counts of [x_n | n]
            for i in range(self.control_q):
                for key in counts:
                    if key[i] == "0":
                        measurement_counts[i][0] += counts[key]
                    else:
                        measurement_counts[i][1] += counts[key]
                                    
            # Recover counts
            for key in counts:
                counts_pp[key] = counts[key]
                for i in range(self.control_q):
                    counts_pp[key] = counts_pp[key] * ( (1 - (eta*self.shots)/(2*measurement_counts[i][int(key[i])])) / (1-eta) )
                counts_pp[key] = round(counts_pp[key])
            counts = counts_pp
        
        return plot_histogram(counts)

    def QpexMod21(self):
        from qiskit import transpile, assemble
        from qiskit.transpiler import Layout
        
        qc = self.ConstructQC() # construct quantum circuit
        backend = self.backend_name  # setting backend        
        t_qc = transpile(qc, backend, optimization_level=2)  # quantum circuit transpile
        result = backend.run(assemble(t_qc, shots=1, memory=1), memory=True).result() # get result
        
        readings = result.get_memory()
        print("Register Reading: " + readings[0])
        
        phase = int(readings[0],2)/(2**self.control_q)
        print("Corresponding Phase: %f" % phase)
        
        return phase

    def RunShorAlgorithm(self):
        from fractions import Fraction
        from math import gcd
            
        factor_found = False
        attempt = 0
        while not factor_found:
            attempt += 1
            print("======================================================")
            
            phase = self.QpexMod21() # Phase = s/r
            print("\nAttempt %i:" % attempt)
            
            frac = Fraction(phase).limit_denominator(self.n) # Denominator should (hopefully!) tell us r
            r = frac.denominator
            print("Result: r = %i" % r)
            
            if phase != 0:
                guesses = [gcd(self.x**(r//2)-1, self.n), gcd(self.x**(r//2)+1, self.n)]
                print("Guessed Factors: %i and %i" % (guesses[0], guesses[1]))
                for guess in guesses:
                    if guess not in [1,self.n] and (self.n % guess) == 0: # Check to see if guess is a factor
                        print("*** Non-trivial factor found: %i ***" % guess)
                        factor_found = True   


# In[3]:


class ShorN35X16:
    
    def __init__(self, n=35, x=16, control_q=3, shots=1000, ancilla_q=2, backend_name="ibmq_qasm_simulator", IsTwirl=False):        
        from math import pi
        
        # load IBMQ Accounts
        backend_use = BackendInit(backend_name)

        # variables
        self.n = n # factoring number
        self.x = x # modular
        self.control_q = control_q # number of control qubits
        self.ancilla_q = ancilla_q # number of ancilla qubits
        self.shots = shots # number of shots
        self.backend_name = backend_use # simulation backend
        self.IsTwirl = IsTwirl # measurement error mitigation (twriling)
        self.Unitary_Twirl = {
            1:[0, 0, 0],
            2:[pi, 3*pi/2, pi/2],
            3:[pi, 0, pi],
            4:[0, pi/2, pi/2],
            5:[pi/2, 0, pi/2],
            6:[pi/2, pi, pi/2],
            7:[pi/2, pi/2, pi],
            8:[pi/2, pi/2, 0],
            9:[pi/2, 3*pi/2, pi],
            10:[pi/2, pi/2, 0],
            11:[pi/2, 0, 3*pi/2],
            12:[pi/2, pi, 3*pi/2]} # twirling parameters
    
    def ConstructQC(self, unitary_index=1):
        from qiskit import QuantumCircuit
        
        qc = QuantumCircuit(self.control_q+self.ancilla_q, self.control_q) # Circuit initialization
        for q in range(self.control_q): # Initialize control qubits
            qc.h(q)

        # Do controlled-U operations
        # Controlled U^1
        qc.cx(2, 4)

        # Controlled U^2
        qc.cx(1, 4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (1, 3, 4))
        qc.cx(4, 3)

        # Controlled U^4
        qc.x(4)
        qc.append(Margolous(self.control_q), (0, 4, 3))
        qc.x(4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (0, 3, 4))
        qc.cx(4, 3)

        # Do inverse-QFT
        qc.append(QftDagger(self.control_q), range(self.control_q))

        if self.IsTwirl == True: # twirling
            qc.append(Twirling(self.Unitary_Twirl, unitary_index, self.control_q), [i for i in range(self.control_q)])
        
        # Measure circuit
        qc.measure(range(self.control_q), range(self.control_q))
        
        return qc
    
    def SimulationQC(self, unitary_index_set=[1,5,9], eta=0.1):
        from qiskit import transpile, assemble
        from qiskit.visualization import plot_histogram
        
        if self.IsTwirl == False: # no twirling
            print('==================================================')
            print('Name of Backend :', self.backend_name) # Setting Simulation Backend
            print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number
            
            qc = self.ConstructQC() # construct quantum circuit
            backend = self.backend_name # setting backend
            t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile

            print("Transpiled Quantum Circuit Depth :", t_qc.depth())
            results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result
            
            # Counts Result
            counts = results.get_counts() # get counts
            print('Simulation Result (counts) :', counts)
            print('==================================================')

        else: # twirling
            ## Quantum pre-processing
            # Variable setting
            counts = {}
            
            # Apply 3 unitray indexes to twriling
            for unitary_index in unitary_index_set:
                print('==================================================')
                print('Unitary Index :', unitary_index) # Setting Unitary Index
                print('Name of Backend :', self.backend_name) # Setting Simulation Backend
                print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number

                qc = self.ConstructQC(unitary_index) # construct quantum circuit
                backend = self.backend_name # setting backend
                t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile
                
                print('Transpiled Quantum Circuit Depth :', t_qc.depth())
                results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result

                # Counts Result
                temp_counts = results.get_counts()
                print('Simulation Result (counts) :', temp_counts)
                
                counts = SumDictionary(counts, temp_counts)
                print('==================================================\n')

            ## Classical post-processing
            # Variable setting
            measurement_counts = []
            for i in range(self.control_q):
                measurement_counts.append([0, 0])
            counts_pp = {}
            
            for key in counts: # average count results
                counts[key] = round(counts[key]/3)
                
            # Calculate counts of [x_n | n]
            for i in range(self.control_q):
                for key in counts:
                    if key[i] == "0":
                        measurement_counts[i][0] += counts[key]
                    else:
                        measurement_counts[i][1] += counts[key]
                                    
            # Recover counts
            for key in counts:
                counts_pp[key] = counts[key]
                for i in range(self.control_q):
                    counts_pp[key] = counts_pp[key] * ( (1 - (eta*self.shots)/(2*measurement_counts[i][int(key[i])])) / (1-eta) )
                counts_pp[key] = round(counts_pp[key])
            counts = counts_pp
        
        return plot_histogram(counts)

    def QpexMod21(self):
        from qiskit import transpile, assemble
        from qiskit.transpiler import Layout
        
        qc = self.ConstructQC() # construct quantum circuit
        backend = self.backend_name  # setting backend        
        t_qc = transpile(qc, backend, optimization_level=2)  # quantum circuit transpile
        result = backend.run(assemble(t_qc, shots=1, memory=1), memory=True).result() # get result
        
        readings = result.get_memory()
        print("Register Reading: " + readings[0])
        
        phase = int(readings[0],2)/(2**self.control_q)
        print("Corresponding Phase: %f" % phase)
        
        return phase

    def RunShorAlgorithm(self):
        from fractions import Fraction
        from math import gcd
            
        factor_found = False
        attempt = 0
        while not factor_found:
            attempt += 1
            print("======================================================")
            
            phase = self.QpexMod21() # Phase = s/r
            print("\nAttempt %i:" % attempt)
            
            frac = Fraction(phase).limit_denominator(self.n) # Denominator should (hopefully!) tell us r
            r = frac.denominator
            print("Result: r = %i" % r)
            
            if phase != 0:
                guesses = [gcd(self.x**(r//2)-1, self.n), gcd(self.x**(r//2)+1, self.n)]
                print("Guessed Factors: %i and %i" % (guesses[0], guesses[1]))
                for guess in guesses:
                    if guess not in [1,self.n] and (self.n % guess) == 0: # Check to see if guess is a factor
                        print("*** Non-trivial factor found: %i ***" % guess)
                        factor_found = True   


# In[4]:


class ShorN65X16:
    
    def __init__(self, n=65, x=16, control_q=3, shots=1000, ancilla_q=2, backend_name="ibmq_qasm_simulator", IsTwirl=False):        
        from math import pi
        
        # load IBMQ Accounts
        backend_use = BackendInit(backend_name)

        # variables
        self.n = n # factoring number
        self.x = x # modular
        self.control_q = control_q # number of control qubits
        self.ancilla_q = ancilla_q # number of ancilla qubits
        self.shots = shots # number of shots
        self.backend_name = backend_use # simulation backend
        self.IsTwirl = IsTwirl # measurement error mitigation (twriling)
        self.Unitary_Twirl = {
            1:[0, 0, 0],
            2:[pi, 3*pi/2, pi/2],
            3:[pi, 0, pi],
            4:[0, pi/2, pi/2],
            5:[pi/2, 0, pi/2],
            6:[pi/2, pi, pi/2],
            7:[pi/2, pi/2, pi],
            8:[pi/2, pi/2, 0],
            9:[pi/2, 3*pi/2, pi],
            10:[pi/2, pi/2, 0],
            11:[pi/2, 0, 3*pi/2],
            12:[pi/2, pi, 3*pi/2]} # twirling parameters
    
    def ConstructQC(self, unitary_index=1):
        from qiskit import QuantumCircuit
        
        qc = QuantumCircuit(self.control_q+self.ancilla_q, self.control_q) # Circuit initialization
        for q in range(self.control_q): # Initialize control qubits
            qc.h(q)

        # Do controlled-U operations
        # Controlled U^1
        qc.cx(2, 4)

        # Controlled U^2
        qc.cx(1, 4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (1, 3, 4))
        qc.cx(4, 3)

        # Controlled U^4
        qc.x(4)
        qc.append(Margolous(self.control_q), (0, 4, 3))
        qc.x(4)
        qc.cx(4, 3)
        qc.append(Margolous(self.control_q), (0, 3, 4))
        qc.cx(4, 3)

        # Do inverse-QFT
        qc.append(QftDagger(self.control_q), range(self.control_q))

        if self.IsTwirl == True: # twirling
            qc.append(Twirling(self.Unitary_Twirl, unitary_index, self.control_q), [i for i in range(self.control_q)])
        
        # Measure circuit
        qc.measure(range(self.control_q), range(self.control_q))
        
        return qc
    
    def SimulationQC(self, unitary_index_set=[1,5,9], eta=0.1):
        from qiskit import transpile, assemble
        from qiskit.visualization import plot_histogram
        
        if self.IsTwirl == False: # no twirling
            print('==================================================')
            print('Name of Backend :', self.backend_name) # Setting Simulation Backend
            print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number
            
            qc = self.ConstructQC() # construct quantum circuit
            backend = self.backend_name # setting backend
            t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile

            print("Transpiled Quantum Circuit Depth :", t_qc.depth())
            results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result
            
            # Counts Result
            counts = results.get_counts() # get counts
            print('Simulation Result (counts) :', counts)
            print('==================================================')

        else: # twirling
            ## Quantum pre-processing
            # Variable setting
            counts = {}
            
            # Apply 3 unitray indexes to twriling
            for unitary_index in unitary_index_set:
                print('==================================================')
                print('Unitary Index :', unitary_index) # Setting Unitary Index
                print('Name of Backend :', self.backend_name) # Setting Simulation Backend
                print('Number of Total Shots :', self.shots) # Setting Simulation Shot Number

                qc = self.ConstructQC(unitary_index) # construct quantum circuit
                backend = self.backend_name # setting backend
                t_qc = transpile(qc, self.backend_name, optimization_level=2) # quantum circuit transpile
                
                print('Transpiled Quantum Circuit Depth :', t_qc.depth())
                results = backend.run(assemble(t_qc, shots=self.shots)).result() # get result

                # Counts Result
                temp_counts = results.get_counts()
                print('Simulation Result (counts) :', temp_counts)
                
                counts = SumDictionary(counts, temp_counts)
                print('==================================================\n')

            ## Classical post-processing
            # Variable setting
            measurement_counts = []
            for i in range(self.control_q):
                measurement_counts.append([0, 0])
            counts_pp = {}
            
            for key in counts: # average count results
                counts[key] = round(counts[key]/3)
                
            # Calculate counts of [x_n | n]
            for i in range(self.control_q):
                for key in counts:
                    if key[i] == "0":
                        measurement_counts[i][0] += counts[key]
                    else:
                        measurement_counts[i][1] += counts[key]
                                    
            # Recover counts
            for key in counts:
                counts_pp[key] = counts[key]
                for i in range(self.control_q):
                    counts_pp[key] = counts_pp[key] * ( (1 - (eta*self.shots)/(2*measurement_counts[i][int(key[i])])) / (1-eta) )
                counts_pp[key] = round(counts_pp[key])
            counts = counts_pp
        
        return plot_histogram(counts)

    def QpexMod21(self):
        from qiskit import transpile, assemble
        from qiskit.transpiler import Layout
        
        qc = self.ConstructQC() # construct quantum circuit
        backend = self.backend_name  # setting backend        
        t_qc = transpile(qc, backend, optimization_level=2)  # quantum circuit transpile
        result = backend.run(assemble(t_qc, shots=1, memory=1), memory=True).result() # get result
        
        readings = result.get_memory()
        print("Register Reading: " + readings[0])
        
        phase = int(readings[0],2)/(2**self.control_q)
        print("Corresponding Phase: %f" % phase)
        
        return phase

    def RunShorAlgorithm(self):
        from fractions import Fraction
        from math import gcd
            
        factor_found = False
        attempt = 0
        while not factor_found:
            attempt += 1
            print("======================================================")
            
            phase = self.QpexMod21() # Phase = s/r
            print("\nAttempt %i:" % attempt)
            
            frac = Fraction(phase).limit_denominator(self.n) # Denominator should (hopefully!) tell us r
            r = frac.denominator
            print("Result: r = %i" % r)
            
            if phase != 0:
                guesses = [gcd(self.x**(r//2)-1, self.n), gcd(self.x**(r//2)+1, self.n)]
                print("Guessed Factors: %i and %i" % (guesses[0], guesses[1]))
                for guess in guesses:
                    if guess not in [1,self.n] and (self.n % guess) == 0: # Check to see if guess is a factor
                        print("*** Non-trivial factor found: %i ***" % guess)
                        factor_found = True   


# In[5]:


def BackendInit(backend_input):
    from qiskit import IBMQ
    
    # backend setting
    provider = IBMQ.load_account()
    my_provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')
    #my_provider = IBMQ.get_provider(hub='ibm-q-skku', group='kaist', project='kaist-graduate')
    backends = my_provider.backends()

    # Check the availability of backend_name
    error_count = 0
    for backend in backends:
        if backend_input == backend.name():
            backend_use = my_provider.get_backend(backend_input)
            break
        else:
            error_count += 1
            
    # Alert Errors...
    if error_count == len(backends):
        raise ValueError("Backend_name must be in the list of IBMQ account backends...")
        
    return backend_use


# In[6]:


BackendInit("ibmq_qasm_simulator")


# In[6]:


def CxMod6(x, power):
    from qiskit import QuantumCircuit
    
    """This function is only valid for x==5"""
    if x not in [5]:
        raise ValueError("'x' must be 5")
    U = QuantumCircuit(3)        

    for iteration in range(power):
        "ccnot-Gate"
        U.x(1)
        U.ccx(1, 2, 0)
        U.x(1)

        "cSWAP-Gate"
        U.x(2)
        U.cswap(2, 0, 1)
        U.x(2)

    U = U.to_gate()
    U.name = "%i^%i mod 6" % (x, power)
    c_U = U.control()
    return c_U    


# In[7]:


def QftDagger(n):
    from qiskit import QuantumCircuit
    import numpy as np
    
    """n-qubit QFTdagger the first n qubits in circ"""
    qc = QuantumCircuit(n)
    # Don't forget the Swaps!
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-np.pi/float(2**(j-m)), m, j)
        qc.h(j)
    qc.name = "QFT†"
    return qc


# In[8]:


def Margolous(n):
    from qiskit import QuantumCircuit
    import numpy as np
    
    """Margolous gate in circ (n=3)"""
    if n != 3:
        raise ValueError("Qubit number of Margolous gate must be n=3.")
    qc = QuantumCircuit(n)
    
    qc.ry(np.pi/4, 2)
    qc.cx(1, 2)
    qc.ry(np.pi/4, 2)
    qc.cx(0, 2)
    qc.ry(-np.pi/4, 2)
    qc.cx(1, 2)
    qc.ry(-np.pi/4, 2)
    
    qc.name = "Marg"
    return qc


# In[1]:


def Twirling(Unitary_Twirl, unitary_index, qubit_num):
    from qiskit import QuantumCircuit, QuantumRegister

    # Create state initialize circuit
    qr = QuantumRegister(qubit_num)
    circ = QuantumCircuit(qr, name='Twirling')
    
    # Apply U
    theta = (Unitary_Twirl[unitary_index])[0]
    phi = (Unitary_Twirl[unitary_index])[1]
    lam = (Unitary_Twirl[unitary_index])[2]
    circ.u(theta, phi, lam, qr)
    circ.barrier()
    
    # Apply U^{dag}
    theta = -(Unitary_Twirl[unitary_index])[0]
    phi = -(Unitary_Twirl[unitary_index])[2]
    lam = -(Unitary_Twirl[unitary_index])[1]
    circ.u(theta, phi, lam, qr)
    
    # Convert to a gate and stick it into an arbitrary place in the bigger circuit
    inst = circ.to_instruction()

    return inst


# In[10]:


def SumDictionary(totdic, adddic):
    for key in adddic:
        if key in totdic:
            totdic[key] += adddic[key]
        else:
            totdic[key] = adddic[key]
    return totdic

