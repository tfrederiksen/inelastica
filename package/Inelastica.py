"""
Inelastica:
Calculates the Inelastic electron tunneling spectra IETS from e-ph coupling
and elastic transport properties (TranSIESTA). The calculations are performed
in the LOE approximation, Paulsson et al., PRB 72, 201101R (2005).
"""
import pyTBT
import SiestaIO as SIO
import EigenChannels as EC
import MiscMath as MM
import numpy as N
import Scientific.IO.NetCDF as NC
import time

mm = MM.mm
dagger = MM.dagger
fermi = MM.fermi
box = MM.box
interpolate = MM.interpolate

class parameter: pass

def main():
    EC.setupParameters(calledFromInelastica=True)

    # Copy calculation parameters ... 
    parameter.Temp, parameter.biasPoints = EC.general.Temp, EC.general.biasPoints 
    parameter.minBias, parameter.maxBias = EC.general.minBias, EC.general.maxBias 
    parameter.modeCutoff, parameter.Ef   = EC.general.modeCutoff, EC.general.Ef 
    parameter.PhHeating, parameter.Vrms  = EC.general.PhHeating, EC.general.Vrms 
    parameter.PhExtDamp                  = EC.general.PhExtDamp 

    parameter.includeModes  = True            # 'True' or range of indices TODO!
    parameter.kT            = parameter.Temp/11604.0 # (eV) 

    EC.readHS()
    readPhononEnergies()
    hw, T, nHT, HT, lLOE, nPhtot, nPh, gamma_eh, gamma_heat = calcIETS()

    V, I, dI, ddI, BdI, BddI, NnPhtot,NnPh = Broaden(lLOE[0],lLOE[1]+lLOE[2],nPhtot,nPh)

    EC.general.iChan, EC.general.iSide = 0, 2
    datafile=EC.fileName()+'.IN'
    initncfile(datafile,hw)
    writeLOEData2Datafile(datafile,hw,T,nHT,HT) 
    writeLOE2ncfile(datafile,hw,nHT,HT,V,I,NnPhtot,NnPh,\
                        dI,ddI,BdI,BddI,gamma_eh,gamma_heat)

def readPhononEnergies():
    NCfile = NC.NetCDFFile(EC.general.PhononNetCDF,'r')
    parameter.hw = N.array(NCfile.variables['hw'][:])
    NCfile.close()


def calcIETS():
    # Calculate the IETS using much of the eigenchannel script for the traces
    EC.calcT(calledFromInelastica=True,InelasticaEf=parameter.Ef)

    # Set up grid and Hilbert term
    kT = parameter.kT 
    
    # Read Self-energies
    T=EC.IETS.T[0][0]

    Egrid=genEgrid() # Energy grid
    hilb=[] # List of hilbert transforms
    ker=None
    Omega=parameter.hw

    # Calculate the prefactors for the Hilbert and non-Hilbert terms
    # Also calculate the Hilbert transfer of the box on the grid for each mode
    for ii in range(len(Omega)):        
        hw=Omega[ii]
        tmp=box(hw,-hw,Egrid,kT)
        tmp2, ker = MM.Hilbert(tmp,ker)
        hilb.append(tmp2)
        SIO.printDone(ii,len(Omega),'LOE : Hilbert transform')

    NN = parameter.biasPoints 

    # Add some points for the Lock in broadening
    approxdV=(parameter.maxBias-parameter.minBias)/NN
    NN+=int(((8*parameter.Vrms)/approxdV)+.5)

    Vl=parameter.minBias-4*parameter.Vrms+ \
        (parameter.maxBias-parameter.minBias+8*parameter.Vrms)/NN*N.array(range(NN),N.float) 

    InH=N.zeros((NN,),N.complex) # Non-Hilb-term current 
    IH=N.zeros((NN,),N.complex) # Hilb-term current
    Pow=N.zeros((NN,),N.float) # Power
    nPhtot=N.zeros((NN,),N.float) # Number of phonons (sum over modes)
    nPhvsBias=N.zeros((NN,len(Omega)),N.float) # Number of phonons 
    nPh=N.zeros((len(Omega),),N.float)         # Number of phonons 
   
    for iV in range(len(Vl)):
        SIO.printDone(iV,len(Vl),'LOE : Calc current')

        InH[iV], IH[iV], Pow[iV] = 0.0, 0.0, 0.0
        V=Vl[iV]

        kasse=box(0,-V,Egrid,kT)
        
        for iOmega in range(len(Omega)):
            hw = Omega[iOmega]
            if hw>parameter.modeCutoff:
                PV=abs(V)

                # Power
                exphw=N.exp(N.clip(hw/kT,-300,300))
                expV=N.exp(N.clip(PV/kT,-300,300))
                exphwmV=N.exp(N.clip((hw-PV)/kT,-300,300))
                if hw/kT>300:
                    mexphw=0.0
                else:
                    mexphw=N.exp(N.clip(-hw/kT,-300,300))
                if PV/kT>300:
                    mexpV=0.0
                else:
                    mexpV=N.exp(N.clip(-PV/kT,-300,300))
                
                # Re-written third time lucky? Numerical problems galore! 
                if (hw-PV)/kT>150:
                    tmpheat=0
                else:
                    if abs((hw-PV)/kT)<1e-3:
                        tmpheat= (2*mexphw*hw+kT*mexphw**2-kT)/(mexphw**2-1)-\
                            (PV-hw)*(4*hw*mexphw**2+kT*mexphw**4-kT)/(2*kT*(1-mexphw**2)**2)
                    else:
                        tmpheat=(1-mexpV)/(1-mexphw)* \
                            ((1+mexphw)*(1-mexpV)*hw-(1-mexphw)*(1+mexpV)*PV)/ \
                            (exphwmV+mexphw*mexpV-mexpV**2-1)
                            
                heat=tmpheat*EC.IETS.P2T[iOmega]/N.pi # Heating term / (hbar hw)

                if parameter.PhHeating: # Heating?
                    nPh[iOmega]=heat/(hw*EC.IETS.P1T[iOmega]/N.pi+parameter.PhExtDamp)+1.0/(exphw-1.0)
                nPhvsBias[iV,iOmega]=nPh[iOmega]

                # Damping term /(hbar hw)
                damp = hw*(1/(exphw-1)-nPh[iOmega])*EC.IETS.P1T[iOmega]/N.pi
            
                # Power in units of 1/(hbar)
                Pow[iV]=hw*(heat+damp)
                nPhtot[iV]=nPhtot[iV]+nPh[iOmega]
                tmp=0.0
                if abs(V-hw)/kT<1e-7:
                    tmp=-kT
                else:
                    tmp=(V-hw)/(N.exp(N.clip((hw-V)/kT,-70.0,70.0))-1)
                if abs(V+hw)/kT<1e-7:
                    tmp+=kT
                else:
                    tmp+=(V+hw)/(N.exp(N.clip((hw+V)/kT,-70.0,70.0))-1)
                InH[iV]+=(tmp-V*nPh[iOmega])*EC.IETS.nHT[iOmega]
                
                # Finite temp Hilbert
                IH[iV]-=EC.IETS.HT[iOmega]*MM.trapez(Egrid,kasse*hilb[iOmega],\
                                                      equidistant=True)/2

        InH[iV]+=V*T # Last to reduce cancelation errors

    # Get the right units for gamma_eh, gamma_heat
    hbar=1.0545727e-34/1.6021773e-19 # hbar in eV S
    gamma_eh=N.zeros((len(Omega),),N.float)
    gamma_heat=N.zeros((len(Omega),),N.float)
    for iOmega in range(len(Omega)):
        # Units [Phonons per Second per dN where dN is number extra phonons]
        gamma_eh[iOmega]=EC.IETS.P1T[iOmega]*Omega[iOmega]/N.pi/hbar
        # Units [Phonons per second per eV [eV-hw]
        gamma_heat[iOmega]=EC.IETS.P2T[iOmega]/N.pi/hbar

    return Omega, T, EC.IETS.nHT, EC.IETS.HT, [Vl, InH, IH], nPhtot, \
           nPhvsBias, gamma_eh, gamma_heat
    

def genEgrid():
    # Generate grid for numerical integration of Hilbert term

    kT=parameter.kT 
    max_hw=max(parameter.hw)
    max_win=max(-parameter.minBias,max_hw)+20*kT+4*parameter.Vrms
    min_win=min(-parameter.maxBias,-max_hw)-20*kT-4*parameter.Vrms

    pts=int(N.floor((max_win-min_win)/kT*3)) 
    grid=N.array(range(pts),N.float)/pts*(max_win-min_win)+min_win

    print "LOE: Hilbert integration grid : %i pts [%f,%f]" % (pts,min(grid),max(grid))

    return grid

################################################################
# Broadening due to Vrms
################################################################

def Broaden(VV,II,nPhtot,nPh):
    """
    Broadening corresponding to Lock in measurements for the
    conductance and IETS spectra. Also resample II, nPh and nPhtot
    to match a common voltage list
    """

    II=II.copy()
    II=II.real
    
    # First derivative dI and bias list dV
    dI=(II[1:len(II)]-II[:-1])/(VV[1]-VV[0])
    dV=(VV[1:len(VV)]+VV[:-1])/2
    # Second derivative and bias ddV
    ddI=(dI[1:len(dI)]-dI[:-1])/(VV[1]-VV[0])
    ddV=(dV[1:len(dV)]+dV[:-1])/2

    # Modulation amplitude
    VA=N.sqrt(2.0)*parameter.Vrms 

    # New bias ranges for broadening
    tmp=int(N.floor(VA/(dV[1]-dV[0]))+1)
    BdV=dV[tmp:-tmp]
    BddV=ddV[tmp:-tmp]

    # Initiate derivatives
    BdI=0*BdV
    BddI=0*BddV
    
    # Calculate first derivative with Vrms broadening
    for iV, V in enumerate(BdV):
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        dIL=interpolate(VL,dV,dI)
        BdI[iV]=2/N.pi*N.sum(dIL*(N.cos(wt)**2))*(wt[1]-wt[0])

    # Calculate second derivative with Vrms broadening    
    for iV, V in enumerate(BddV):
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        ddIL=interpolate(VL,ddV,ddI)
        BddI[iV]=8.0/3.0/N.pi*N.sum(ddIL*(N.cos(wt)**4))*(wt[1]-wt[0])

    # Reduce to one voltage grid
    NN=parameter.biasPoints 
    V=parameter.minBias+(parameter.maxBias-parameter.minBias)/NN*N.array(range(NN)) 

    NI=interpolate(V,VV,II)
    NdI=interpolate(V,dV,dI)
    NddI=interpolate(V,ddV,ddI)
    NBdI=interpolate(V,BdV,BdI)
    NBddI=interpolate(V,BddV,BddI)
    NnPhtot=interpolate(V,VV,nPhtot)
    NnPh=N.zeros((len(V),len(nPh[0,:])),N.float)
    for ii in range(len(nPh[0,:])):
        NnPh[:,ii]=interpolate(V,VV,nPh[:,ii])
    
    return V, NI ,NdI, NddI, NBdI, NBddI, NnPhtot, NnPh

################################################################
# Output to NetCDF file
################################################################

def initncfile(filename,hw):
    'Initiate netCDF file'
    print 'Inelastica: Initializing nc-file'
    ncfile = NC.NetCDFFile(filename+'.nc','w','Created '+time.ctime(time.time()))
    ncfile.title = 'Inelastica Output'
    ncfile.version = 1
    ncfile.createDimension('Nph',len(hw))
    ncfile.createDimension('Bias',None)
    tmp = ncfile.createVariable('IelasticL','d',('Bias',))
    tmp.units = 'Conductance quantum'
    tmp1 = ncfile.createVariable('IelasticR','d',('Bias',))
    tmp1.units = 'Conductance quantum'
    tmp2 = ncfile.createVariable('IinelasticL','d',('Bias',))
    tmp2.units = 'Conductance quantum'
    tmp3 = ncfile.createVariable('IinelasticR','d',('Bias',))
    tmp3.units = 'Conductance quantum'
    tmp4 = ncfile.createVariable('Bias','d',('Bias',))
    tmp4.units = 'V'
    nchw = ncfile.createVariable('hw','d',('Nph',))
    nchw[:] = N.array(hw)
    nchw.units = 'eV'
    ncfile.close()
    
def writeLOE2ncfile(filename,hw,nHT,HT,V,I,nPhtot,nPh,\
                     dI,ddI,BdI,BddI,gamma_eh,gamma_heat):
    'Write LOE data to netCDF file'
    print 'Inelastica: Write LOE data to nc-file'
    ncfile = NC.NetCDFFile(filename+'.nc','a')
    ncfile.createDimension('LOE_V',len(V))
    tmp=ncfile.createVariable('LOE_V','d',('LOE_V',))
    tmp[:]=N.array(V)
    tmp.units='V'
    tmp=ncfile.createVariable('LOE_IETS','d',('LOE_V',))
    tmpn=BddI/BdI
    tmp[:]=N.array(tmpn.real)
    tmp.units='Broadened ddI/dI [1/V]'
    tmp=ncfile.createVariable('LOE_dI','d',('LOE_V',))
    tmp[:]=N.array(dI.real)
    tmp.units='dI LOE, transmission/V'
    tmp=ncfile.createVariable('LOE_ddI','d',('LOE_V',))
    tmp[:]=N.array(ddI.real)
    tmp.units='ddI LOE, transmission/V^2'
    tmp=ncfile.createVariable('LOE_BdI','d',('LOE_V',))
    tmp[:]=N.array(BdI.real)
    tmp.units='Broadened dI LOE, transmission/V'
    tmp=ncfile.createVariable('LOE_BddI','d',('LOE_V',))
    tmp[:]=N.array(BddI.real)
    tmp.units='Broadened ddI LOE, transmission/V^2'
    tmp=ncfile.createVariable('LOE_I','d',('LOE_V',))
    tmp[:]=N.array(I.real)
    tmp.units='transmission/V'
    tmp=ncfile.createVariable('LOE_tot_nPh','d',('LOE_V',))
    tmp[:]=N.array(nPhtot.real)
    tmp.units='Total number of phonons '
    tmp=ncfile.createVariable('LOE_nPh','d',('LOE_V','Nph'))
    tmp[:]=N.array(nPh.real)
    tmp.units='Number of phonons'
    tmp=ncfile.createVariable('LOE_gamma_eh','d',('Nph',))
    tmp[:]=N.array(gamma_eh)
    tmp.units='e-h damping [*deltaN=1/Second]'
    tmp=ncfile.createVariable('LOE_gamma_heat','d',('Nph',))
    tmp[:]=N.array(gamma_heat)
    tmp.units='Phonon heating [*(bias-hw) (eV) = 1/Second]'
    unitConv=1.602177e-19/N.pi/1.054572e-34
    ncTotR = ncfile.createVariable('LOE_ISymTr','d',('Nph',))
    ncTotR[:] = N.array(nHT)*unitConv
    ncTotR.units = 'Trace giving Symmetric current contribution (same units as eigenchannels [1/s/eV])'
    ncTotL = ncfile.createVariable('LOE_IAsymTr','d',('Nph',))
    ncTotL[:] = N.array(HT)
    ncTotL.units = 'Trace giving Asymmetric current contribution'
    
    ncfile.close()
    
def writeLOEData2Datafile(file,hw,T,nHT,HT):
    f = open(file,'a')
    f.write("## Almost First Order Born Approximation (kbT=0)\n")
    f.write("## hw(eV)      Trans        non-H-term (e/pi/hbar)   H-term\n")
    for ii in range(len(hw)):
        f.write("## %e %e %e %e\n" % (hw[ii],T,nHT[ii],HT[ii]))
    f.write("\n")
    f.close()

def writeData2ncfile(filename,dict,ibias):
    'Append data to netCDF file'
    print 'Inelastica: Writing SCBA data to nc-file'
    bias = dict['bias']
    IelasticL = dict['NCL'] # Elastic current [Units, e/Pi/hbar, i.e.,
    #                         energy in eV, bias in V -> Current * 77 muA]
    IelasticR = (-dict['NCR'])
    IinelasticL = (dict['ICL']) # Left Current
    IinelasticR = (-dict['ICR']) # Right Current
    Nph=dict['nPh'] # Number of phonons
    # Add new data (some nicer way should exist)
    
    ncfile = NC.NetCDFFile(filename+'.nc','a')
    ncIelasticL = ncfile.variables['IelasticL']
    ncIelasticR = ncfile.variables['IelasticR']
    ncIinelasticL = ncfile.variables['IinelasticL']
    ncIinelasticR = ncfile.variables['IinelasticR']
    ncBias = ncfile.variables['Bias']
    ncNph = ncfile.variables['SCBA_nPh']

    tmp = ncBias[:]
    ncBias[len(tmp)] = bias

    tmp = ncIelasticL[:]
    ncIelasticL[len(tmp)]=IelasticL
    ncIelasticR[len(tmp)]=IelasticR
    ncIinelasticL[len(tmp)]=IinelasticL
    ncIinelasticR[len(tmp)]=IinelasticR
    ncNph[len(tmp),:]=Nph
    ncfile.close()



def writeData2Datafile(file,dict,hw):
    c1 = dict['bias']
    c2 = dict['NCL'] # Currents in units of e/pi/hbar*energy, i.e.,
    #                  G=77 muA/V * dI/dV if energy [eV], bias [V]
    c3 = (dict['ICL']-dict['ICR'])/2 # Average of left and right
    c4 = dict['NPL']+dict['NPR']
    c5 = dict['IPL']+dict['IPR']
    data = [c1,c2,c3,c4,c5]
    form = '%e   %e   %e   %e   %e   '
    for i in range(len(hw)):
        data.append(dict['nPh'][i])
        form += '%f   '
    form += '\n'
    data = tuple(data)
    f = open(file,'a')
    f.write(form %data)
    f.close()

##################### Start main routine #####################

if __name__ == '__main__':
    main()




