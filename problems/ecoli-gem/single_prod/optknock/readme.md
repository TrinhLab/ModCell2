# Description
Round 1: All models where solved with cplex with a limit of 100 min
Round 2: Those models which could not be solved on round 1 were solved in cplex with a limit of 200 min.
- After round 2 some models failed to solve correctly (Cplex reported an optimal solution, but upon further inspection, e.g. porduction envelope, the solutions would abolish growth and not behanve as expected)
Models where then solved with gurobi, which could only solve 2 correctly. 
Round 3 (file deleted): All models solved with cplex with a time limit of 10,000 seconds (166 min), however this time the default big-M bounds for reactions were left (Original these bounds where tightened through FBA). While the models solved a lot faster, some of them stilled failed (less than half). 
Round 4: Original model bounds and modified integrality tolerance, which was lowered from 1e-5 (default) to 1e-8. 4 knockouts.
