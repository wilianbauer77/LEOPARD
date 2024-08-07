using Pkg
Pkg.add("Plots")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
using CSV
using DataFrames
using Plots
using DelimitedFiles

using DelimitedFiles
using Plots

# Carregar os dados dos arquivos .dat
data1 = readdlm("leopardv3001.dat")
data2 = readdlm("leopardv3002.dat")
data3 = readdlm("leopardv3003.dat")
data4 = readdlm("leopardv3004.dat")
data5 = readdlm("leopardv3005.dat")
data6 = readdlm("leopardv3006.dat")
data7 = readdlm("leopardv3007.dat")
data8 = readdlm("leopardv3008.dat")
data9 = readdlm("leopardv3009.dat")
data10 = readdlm("leopardv3010.dat")

# Extrair as colunas desejadas
x1 = data1[:, 1]
y1 = data1[:, 2]

x2 = data2[:, 1]
y2 = data2[:, 2]

# Criar o gráfico
plot(x1, y1, label="Dados 1", title="Gráfico de Dados Múltiplos", xlabel="Eixo X", ylabel="Eixo Y", marker=:circle)
plot!(x2, y2, label="Dados 2", marker=:square)

# Salvar o gráfico em um arquivo
savefig("grafico_multiplos_dados.png")
