compilando:
    cd ./base_heuristic
    make
    depois
    cd ./proj_heuristic
    make
agora os executaveis foram compilados e podemos rodar run_experiments.sh

o run_experiments gera uma saida "results.txt" que executa tanto a heuristica base 
como a nossa heuristica para todas os dados que estão em "./data".
Em run_experiments tem um loop que varia o parametro alfa (densidade minima das cliques)
e o tempo limite de execução. Ainda, como as duas heuristicas possuim aleatoriedade, tem
o campo seed que controla isso.