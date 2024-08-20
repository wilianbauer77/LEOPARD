import numpy
import mpmath
import matplotlib.pyplot as plt

# Número de espécies
Nspecies = 2

# Valores de kappa para comparação
k = [50, 10, 6, 4]

npara = numpy.zeros(Nspecies, dtype='i4')
nperp = numpy.zeros(Nspecies, dtype='i4')

vparamin = numpy.zeros(Nspecies)
vparamax = numpy.zeros(Nspecies)

vperpmin = numpy.zeros(Nspecies)
vperpmax = numpy.zeros(Nspecies)

dens = numpy.zeros(Nspecies)
mu = numpy.zeros(Nspecies)

beta_para = numpy.zeros(Nspecies)
beta_perp = numpy.zeros(Nspecies)

theta_para = numpy.zeros(Nspecies)
theta_perp = numpy.zeros(Nspecies)
vdrift = numpy.zeros(Nspecies)

# Configurações das espécies
npara[0] = 255
nperp[0] = 64
dens[0] = 1.0
mu[0] = 1.0
beta_para[0] = 1.0
beta_perp[0] = 1.0
theta_para[0] = (((k[0] - 1.5) / k[0]) * beta_para[0] * mu[0] / dens[0]) ** 0.5
theta_perp[0] = (((k[0] - 1.5) / k[0]) * beta_perp[0] * mu[0] / dens[0]) ** 0.5
vdrift[0] = -0.0
vparamin[0] = -5.0
vparamax[0] = 5.0
vperpmin[0] = -10.0
vperpmax[0] = 10.0

npara[1] = 255
nperp[1] = 64
dens[1] = 1.0
mu[1] = 1836.0
beta_para[1] = 20.0
beta_perp[1] = 1.0
theta_para[1] = (((k[1] - 1.5) / k[1]) * beta_para[1] * mu[1] / dens[1]) ** 0.5
theta_perp[1] = (((k[1] - 1.5) / k[1]) * beta_perp[1] * mu[1] / dens[1]) ** 0.5
vdrift[1] = -0.0
vparamin[1] = -300.0
vparamax[1] = 300.0
vperpmin[1] = -130.0
vperpmax[1] = 130.0

limit = 10.0 ** (-300)

def dist_bimax(vpar, vper, n, mu, beta_par, beta_per, drift):
    bimax = numpy.exp(-n * (vpar - drift) ** 2 / beta_par / mu - n * vper ** 2 / beta_per / mu) * n ** 1.5 / (mu ** 1.5 * numpy.pi ** 1.5 * beta_per * numpy.sqrt(beta_par))
    return bimax.ravel()

def dist_bikappa(vpar, vper, n, mu, kappa, theta_par, theta_per, drift):
    theta_par_per2 = numpy.outer(theta_par, theta_per ** 2)
    bikappa = (n / (kappa ** 1.5 * numpy.pi ** 1.5 * theta_par_per2)) * (mpmath.gamma(kappa + 1) / (mpmath.gamma(kappa - 0.5))) * (1 + (vpar - drift) ** 2 / (kappa * theta_par ** 2) + vper ** 2 / (kappa * theta_per ** 2)) ** (-kappa - 1)
    return bikappa.ravel()

for ispecies in range(0, Nspecies):
    file_name_bimax = 'distribution_bimax' + str(ispecies + 1) + '.dat'
    file_name_bikappa = 'distribution_bikappa_' + str(k[ispecies]) + '_' + str(ispecies + 1) + '.dat'

    # Calcula a distribuição bi-maxwelliana
    vpara = numpy.linspace(vparamin[ispecies], vparamax[ispecies], npara[ispecies])
    vperp = numpy.linspace(vperpmin[ispecies], vperpmax[ispecies], nperp[ispecies])
    vpara2, vperp2 = numpy.meshgrid(vpara, vperp)

    databimax = dist_bimax(vpara2, vperp2, dens[ispecies], mu[ispecies], beta_para[ispecies], beta_perp[ispecies], vdrift[ispecies])

    # Calcular o valor médio de v_perp
    vperp_mean = numpy.mean(vperp)

    for kappa in k:
        theta_para[ispecies] = (((kappa - 1.5) / kappa) * beta_para[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5
        theta_perp[ispecies] = (((kappa - 1.5) / kappa) * beta_perp[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5

        databikappa = dist_bikappa(vpara2, vperp2, dens[ispecies], mu[ispecies], kappa, theta_para[ispecies], theta_perp[ispecies], vdrift[ispecies])
        databikappa = databikappa.reshape(nperp[ispecies], npara[ispecies])

        data_new_bikappa = numpy.zeros((nperp[ispecies], npara[ispecies]))

        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                if databikappa[j, i] > limit:
                    data_new_bikappa[j, i] = databikappa[j, i] * dens[ispecies]
                else:
                    data_new_bikappa[j, i] = 0.0

        dat_fin_bimax = open(file_name_bimax, 'w')
        dat_fin_bikappa = open(file_name_bikappa, 'w')
        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                dat_fin_bimax.write(str(vpara[i]))
                dat_fin_bikappa.write(str(vpara[i]))
                dat_fin_bimax.write(" ")
                dat_fin_bikappa.write(" ")
                dat_fin_bimax.write(str(vperp[j]))
                dat_fin_bikappa.write(str(vperp[j]))
                dat_fin_bimax.write(" ")
                dat_fin_bikappa.write(" ")
                dat_fin_bimax.write(str(data_new_bikappa[j, i]))
                dat_fin_bikappa.write(str(data_new_bikappa[j, i]))
                dat_fin_bimax.write("\n")
                dat_fin_bikappa.write("\n")

        # Plot da distribuição bi-kappa para o valor atual de kappa
        plt.plot(vpara, numpy.log10(data_new_bikappa[:, :].sum(axis=0)), label=rf'$\kappa={kappa}, v_\perp={vperp_mean:.2f}$')

    plt.xlabel(r'$v_\parallel$')
    plt.ylabel(r'$\log{f(v)}$')
    plt.title(f'Species {ispecies + 1} Distribution')
    plt.legend()
    plt.savefig(file_name_bimax + '.png')  # Salva o gráfico como imagem
    plt.savefig(file_name_bikappa + '.png')  # Salva o gráfico como imagem
    plt.show()