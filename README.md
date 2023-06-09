# multiobjective-reliability-based-design-optimization-algorithm
Multiobjective reliability-based design optimization (RBDO) is a research area, which has not been investigated in the literatures comparing with single-objective RBDO. This work conducts an exhaustive study of fifteen new and popular metaheuristic multiobjective RBDO algorithms, including non-dominated sorting genetic algorithm II, differential evolution for multiobjective optimization, multiobjective evolutionary algorithm based on decomposition, multiobjective particle swarm optimization, multiobjective flower pollination algorithm, multiobjective bat algorithm, multiobjective grey wolf optimizer, multiobjective multi-verse optimization, multiobjective water cycle optimizer, success history-based adaptive multiobjective differential evolution, success history-based adaptive multiobjective differential evolution with whale optimization, multiobjective salp swarm algorithm, real-code population-based incremental learning and differential evolution, unrestricted population size evolutionary multiobjective optimisation algorithm, and multiobjective jellyfish search optimizer. In addition, the adaptive chaos control method is employed for the above-mentioned algorithms to estimate the probabilistic constraints effectively. This comparative analysis reveals the critical technologies and enormous challenges in the RBDO field. It also offers new insight into simultaneously dealing with the multiple conflicting design objectives and probabilistic constraints. Also, this study presents the advantage and future development trends or incurs the increased challenge of researchers to put forward an effective multiobjective RBDO algorithm that assists the complex engineering system design.

The descriptions of the two examples are as follows:
1. Multiobjective RBDO for spring.

The multiobjective RBDO of spiral squeezing spring is selected as the first example, which is modified from (Kannan and Kramer 1994; Kumar et al. 2021a). The multiobjective RBDO of spring can be deemed as a real mechanical application example with integer, discrete, and continuous variables. The schematic view is shown in Fig. 4. There are two design objectives and eight constraints. The first objective is minimizing volume of the steel wire, while the second objective is minimizing the shear stress. The integer variable   and the discrete variable   are considered as the deterministic design variables, while the continuous variable   is selected as the random variable that obeys the normal distribution with standard deviation 0.005.

2. Multiobjective RBDO for simply supported I-beam.

A simply supported I-beam is used as the second examples, which is widely tested by different researches (Huang et al. 2006; Kumar et al. 2021a). The schematic diagram is plotted in Fig. 4. There are two loads, i.e. P=600 kN and Q=50 kN. The Young’s modulus E is 2×104 kN/cm2. The length L is 200 cm. The allowable bending stress σb is 16 kN/cm2. Two design objectives are minimizations of cross-sectional area and vertical deflection. The stress constraint is utilized. The geometric dimensions are selected as the design and random variables, which are assumed following the normal distribution with standard deviation 0.005.

The details of running the program are as follows:
1. Run main.m to obtain the results of these methods.
2. Run main_filter.m to filter the solutions in the infeasible domain.
3. Run main_generate_results_new to obtain the all results.

![image](https://github.com/hfut-mengz/multiobjective-reliability-based-design-optimization-algorithm/blob/main/flowchart.jpg)

For more information, please refer to the following: (https://github.com/hfut-mengz/multiobjective-reliability-based-design-optimization-algorithm.git)
