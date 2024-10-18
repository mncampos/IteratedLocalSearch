using JuMP, Gurobi

function read_instance(filename::String)
    open(filename, "r") do file
        line = readline(file)
        n, m, W = parse.(Int, split(line))
        
        readline(file) 
        
        t = Int[]
        while true
            line = readline(file)
            if line == ""  
                break
            end
            append!(t, parse.(Int, split(line)))
        end
        
        w = Int[]
        while true
            line = readline(file)
            if line == ""
                break
            end
            append!(w, parse.(Int, split(line)))
        end

        I = []
        while !eof(file)
            line = readline(file)
            if line != "" 
                values = parse.(Int, split(line))
                if length(values) == 2
                    push!(I, (values[1], values[2]))
                end
            end
        end
        
        return n, w, t, W, I
    end
end

# Carregar dados
filename = "ep01.dat"  # NOME DO ARQUIVO DA INSTÂNCIA
n, w, t, W, I = read_instance(filename)

# Criar modelo
print("Criando modelo! \n");
m = Model(Gurobi.Optimizer)
set_attribute(m, "TimeLimit", 1800)
set_attribute(m, "Presolve", 2)

print("Criando variáveis e função objetivo! \n");
@variable(m, x[1:n], Bin)
@objective(m, Max, sum(x[i] * t[i] for i in 1:n))

print("Construindo restrições! \n");
@constraint(m, sum(x[i] * w[i] for i in 1:n) <= W)


for (j, k) in I
    @constraint(m, x[j] + x[k] <= 1)
end

print("Iniciando otimização! \n")
optimize!(m)

println("Status: ", termination_status(m))
println("Ingredientes selecionados: ", [i for i in 1:n if value(x[i]) > 0.5])
println("Sabor total: ", objective_value(m))
