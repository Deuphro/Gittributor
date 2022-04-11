#pragma TextEncoding = "UTF-8"
//date: september 1st, 2014
//owner: Francois-Regis ORTHOUS-DAUNAY, working at IPAG (UMR5374), CNRS, UJF in Grenoble, FRANCE
//recipient:

//1 license:
//The owner grants non exclusive and non transferable license to use this software and its 
//documentation only for his or her research activities as well as a copy of ATTRIBUTOR for archive 
//and backup purpose as long as the titles, trademarks and other personal and reserved rights are 
//mentioned on every copies and they are subjected to this license.

//2 distribution:
//the recipient should not:
//1) make usable by anyone or distribute to anyone this software or any element of this software or related documentation.
//2) grant a sub-license or copy, design a reversed version, dissociate or decompose, modify partly or in general this software or its documentation.
//3) use this software in service exchange terms, or subcontracting, in a service environment or make this software accessible for third party.

//3 Property
//The recipient recognizes copyright and trade secret as well as every other rights relative to intellectual property of any kind on this software, its
//documentation and specifications which are all and will be owner's property. The recipient recognizes that none of the items in this license will
//be interpreted as transfering any part of the owneer's right to third party.

#pragma rtGlobals=1		// Use modern global access method.
#include  <SaveGraph>
#include <Resize Controls>
#include <Multi-peak fitting 2.0>

constant emass=0.00054857990946

macro makelesproc()
variable indice=0
string nomproc
do
nomproc="elem"+elements(indice)
print  "Function",nomproc,"(ba) : ButtonControl"
print "	STRUCT WMButtonAction &ba"
print "wave/Z molecule"
print "	switch( ba.eventCode )"
print "		case 2: "
print "			molecule[",num2str(indice),"]+=1"
print "			genestringmol()"
print "			break"
print "		case 3: "
print "			molecule[",num2str(indice),"]-=1"
print "			genestringmol()"
print "			break"
print "	endswitch"
print "	return 0"
print "End"
indice+=1
while (indice<114)
endmacro


Function ButtonReset(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave/Z Molecule
nvar krikriter
	switch( ba.eventCode )
		case 2:
			Molecule=0
			tuelesspec()
			genestringmol()
			labelize(simu_mass,simu_proba)
			//krikriter=0.01
			break
	endswitch

	return 0
End

function genestringmol()
wave/Z molecule, nb_isotopes,mendev_masses, mendev_abundances
wave/T elements
variable n=numpnts(elements), k,i, j, nbmass=0, massmol=0, intmass=0
Nvar krikriter, charge
string elem
string/G nommol
nommol=""
	for(k=0;k<n;k+=1)//verifie la presence de chaque element
		if (molecule[k]!=0) // si l' element k est present
			elem=elements[k]
			nommol=nommol+elem+"\\B"+num2str(molecule[k])+"\M"//genere la formule brute
			choisylecrible(k,molecule[k])
		endif
	endfor
	if (charge>0)
		nommol=nommol+"\S"+num2str(abs(charge))+"+\M"
	elseif (charge<0)
		nommol=nommol+"\S"+num2str(abs(charge))+"-\M"
	endif
	generesimu(krikriter,charge)
	//titlebox title0 variable=nommol, pos={313,1}, labelBack=(65535,65535,65535)
	setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
	setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
	if(stringmatch(nommol[0],"\*")==1)
		ModifyGraph/W=dmvm hideTrace(simu_probadm)=1
		ModifyGraph/W=agregateur hideTrace(simu_proba)=1
	elseif(stringmatch(nommol[0],"\*")==0)
		ModifyGraph/W=dmvm hideTrace(simu_probadm)=0
		ModifyGraph/W=agregateur hideTrace(simu_proba)=0
	endif
end

Function/S mol2str(molecule,charge)
wave molecule
variable charge
wave/T elements
variable n=numpnts(elements), k,t=ticks
string elem
string nommol
nommol=""
	for(k=0;k<n;k+=1)//verifie la presence de chaque element
		if (molecule[k]!=0) // si l' element k est present
			elem=elements[k]
			nommol=nommol+elem+"\\B"+num2str(molecule[k])+"\M"//genere la formule brute
		endif
	endfor
	if (charge>0)
		nommol=nommol+"\S"+num2str(abs(charge))+"+\M"
	elseif (charge<0)
		nommol=nommol+"\S"+num2str(abs(charge))+"-\M"
	endif
return nommol
end

Function/S mol2instr(molecule,charge)
wave molecule
variable charge
wave/T elements
variable n=numpnts(elements), k,t=ticks
string elem
string nommol
nommol=""
	for(k=0;k<n;k+=1)//verifie la presence de chaque element
		if (molecule[k]!=0) // si l' element k est present
			elem=elements[k]
			nommol=nommol+elem+num2str(molecule[k])+" "//genere la formule brute
		endif
	endfor
	if (charge>0)
		nommol+="+"+num2str(charge)
	elseif (charge<0)
		nommol+=num2str(charge)
	endif
return nommol
end

Function mol2massmono(molecule,charge)
wave molecule
variable charge
wave mendev_masses
wave/T elements
variable n=numpnts(elements), k,t=ticks
variable massmol=0
	for(k=0;k<n;k+=1)//verifie la presence de chaque element
		if (molecule[k]!=0) // si l' element k est present
			massmol+=molecule[k]*mendev_masses[k][0]
		endif
	endfor
	massmol=(massmol-charge*emass)/max(abs(charge),1)
return massmol
end

function crible10(e,n)
variable n,e
variable k0,k1,k2,k3,k4,k5,k6,k7,k8,k9, count, counbi=0, i=10,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor
	for(k0=0;k0<n+1;k0+=1)	//pour k0 fois le premier isotope il reste n+1-k0 places disponible
		for(k1=0;k1<n+1-k0;k1+=1)	// pour k1 fois le deuxieme ET k0 fois le premier il reste n-k0-k1 places
			for(k2=0;k2<n+1-k1-k0;k2+=1)
				for(k3=0;k3<n+1-k2-k1-k0;k3+=1)
					for(k4=0;k4<n+1-k3-k2-k1-k0;k4+=1)
						for(k5=0;k5<n+1-k4-k3-k2-k1-k0;k5+=1)	
							for(k6=0;k6<n+1-k5-k4-k3-k2-k1-k0;k6+=1)
								for(k7=0;k7<n+1-k6-k5-k4-k3-k2-k1-k0;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5-k4-k3-k2-k1-k0;k8+=1)	
											k9=n-k8-k7-k6-k5-k4-k3-k2-k1-k0
											combisotopes[0]={k0,k1,k2,k3,k4,k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
					endfor
				endfor
			endfor
		endfor
	endfor
end

function crible9(e,n)
variable n,e
variable k1,k2,k3,k4,k5,k6,k7,k8,k9, count, counbi=0, i=9,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		for(k1=0;k1<n+1;k1+=1)	// pour k1 fois le deuxieme ET k0 fois le premier il reste n-k0-k1 places
			for(k2=0;k2<n+1-k1;k2+=1)
				for(k3=0;k3<n+1-k2-k1;k3+=1)
					for(k4=0;k4<n+1-k3-k2-k1;k4+=1)
						for(k5=0;k5<n+1-k4-k3-k2-k1;k5+=1)	
							for(k6=0;k6<n+1-k5-k4-k3-k2-k1;k6+=1)
								for(k7=0;k7<n+1-k6-k5-k4-k3-k2-k1;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5-k4-k3-k2-k1;k8+=1)	
											k9=n-k8-k7-k6-k5-k4-k3-k2-k1
											combisotopes[0]={k1,k2,k3,k4,k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
					endfor
				endfor
			endfor
		endfor
end

function crible8(e,n)
variable n,e
variable k2,k3,k4,k5,k6,k7,k8,k9, count, counbi=0, i=8,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
			for(k2=0;k2<n+1;k2+=1)
				for(k3=0;k3<n+1-k2;k3+=1)
					for(k4=0;k4<n+1-k3-k2;k4+=1)
						for(k5=0;k5<n+1-k4-k3-k2;k5+=1)	
							for(k6=0;k6<n+1-k5-k4-k3-k2;k6+=1)
								for(k7=0;k7<n+1-k6-k5-k4-k3-k2;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5-k4-k3-k2;k8+=1)	
											k9=n-k8-k7-k6-k5-k4-k3-k2
											combisotopes[0]={k2,k3,k4,k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
					endfor
				endfor
			endfor
end

function crible7(e,n)
variable n,e
variable k3,k4,k5,k6,k7,k8,k9, count, counbi=0, i=7,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
				for(k3=0;k3<n+1;k3+=1)
					for(k4=0;k4<n+1-k3;k4+=1)
						for(k5=0;k5<n+1-k4-k3;k5+=1)	
							for(k6=0;k6<n+1-k5-k4-k3;k6+=1)
								for(k7=0;k7<n+1-k6-k5-k4-k3;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5-k4-k3;k8+=1)	
											k9=n-k8-k7-k6-k5-k4-k3
											combisotopes[0]={k3,k4,k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
					endfor
				endfor
end

function crible6(e,n)
variable n,e
variable k4,k5,k6,k7,k8,k9, count, counbi=0, i=6,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
					for(k4=0;k4<n+1;k4+=1)
						for(k5=0;k5<n+1-k4;k5+=1)	
							for(k6=0;k6<n+1-k5-k4;k6+=1)
								for(k7=0;k7<n+1-k6-k5-k4;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5-k4;k8+=1)	
											k9=n-k8-k7-k6-k5-k4
											combisotopes[0]={k4,k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
					endfor
end

function crible5(e,n)
variable n,e
variable k5,k6,k7,k8,k9, count, counbi=0, i=5,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
						for(k5=0;k5<n+1;k5+=1)	
							for(k6=0;k6<n+1-k5;k6+=1)
								for(k7=0;k7<n+1-k6-k5;k7+=1)
									for(k8=0;k8<n+1-k7-k6-k5;k8+=1)	
											k9=n-k8-k7-k6-k5
											combisotopes[0]={k5,k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
						endfor
end

function crible4(e,n)
variable n,e
variable k6,k7,k8,k9, count, counbi=0, i=4,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
							for(k6=0;k6<n+1;k6+=1)
								for(k7=0;k7<n+1-k6;k7+=1)
									for(k8=0;k8<n+1-k7-k6;k8+=1)	
											k9=n-k8-k7-k6
											combisotopes[0]={k6,k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
							endfor
end

function crible3(e,n)
variable n,e
variable k7,k8,k9, count, counbi=0, i=3,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
								for(k7=0;k7<n+1;k7+=1)
									for(k8=0;k8<n+1-k7;k8+=1)	
											k9=n-k8-k7
											combisotopes[0]={k7,k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
								endfor
end

function crible2(e,n)
variable n,e
variable k8,k9, count, counbi=0, i=2,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

		
									for(k8=0;k8<n+1;k8+=1)	
											k9=n-k8
											combisotopes[0]={k8,k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
									endfor
end

function crible1(e,n)
variable n,e
variable k9, count, counbi=0, i=1,ka, atome= e// chaque isotope a n+1 occurence possible
Killwaves/Z mass, proba,combisotopes, masses, occurence
Make/O/D/N=0 mass, proba
Make/O/D/N=(i) combisotopes, masses, occurence
wave/Z combisotopes, mendev_masses, mendev_abundances, masses, occurence, proba, mass

	for(count=0;count<i;count+=1)
		masses[count]=mendev_masses[atome][count]
		occurence[count]=mendev_abundances[atome][count]
	endfor

											k9=n-k8
											combisotopes[0]={k9}
											InsertPoints counbi,1, proba
											proba[counbi]=1
											InsertPoints counbi,1, mass
											mass[counbi]=0
												for(count=0;count<i;count+=1)
													proba[counbi] = proba[counbi] * binomial(sum(combisotopes,count,i),combisotopes[count]) * occurence[count]^combisotopes[count] //pourchaque isotope il existe 'choix de kx atome parmi le nombre de place restantes'
													mass[counbi] = mass[counbi] + masses[count]*combisotopes[count]
												endfor
											counbi+=1
end

function choisylecrible(element,Zstoechio)
variable element, Zstoechio
variable stoechio=abs(Zstoechio)
string blaz, bapro, ssema
wave/Z nb_isotopes
wave/T elements
blaz=elements[element]
bapro=blaz+"_probabilites"
ssema=blaz+"_massexact"
variable i=nb_isotopes[element]
	if (i==1)
		crible1(element,stoechio)
	elseif (i==2)
		crible2(element,stoechio)
	elseif (i==3)
		crible3(element,stoechio)
	elseif (i==4)
		crible4(element,stoechio)
	elseif (i==5)
		crible5(element,stoechio)
	elseif (i==6)
		crible6(element,stoechio)
	elseif (i==7)
		crible7(element,stoechio)
	elseif (i==8)
		crible8(element,stoechio)
	elseif (i==9)
		crible9(element,stoechio)
	elseif (i==10)
		crible10(element,stoechio)
	endif
wave/Z proba, mass
mass*=sign(Zstoechio)
Duplicate/O proba $bapro
Duplicate/O mass $ssema
end

function tuelesspec()
variable k,n
string lesprobas, lesmasses, cible
lesprobas=Wavelist("*_probabilites",";","")
lesmasses=Wavelist("*_massexact",";","")
n=itemsinlist(lesprobas)
	for(k=0;k<n;k+=1)	// Initialize variables;continue test
		cible=stringfromlist(k,lesprobas)
		Killwaves/Z $cible
		cible=stringfromlist(k,lesmasses)
		Killwaves/Z $cible
	endfor
Killwaves/Z simu_proba, simu_mass, simu_probadm
Make/O/D/N=1 simu_proba, simu_mass, simu_probadm
simu_proba=1
simu_mass=0
simu_probadm=defmass(simu_mass)	
end

function generesimu(kriter,charge) //kriter est le ratio minimum d'une probabilité isotopique pour etre prise en compte dans la formule entre 0 pour tout prendre et 1 pour avoir le plus probable seulement
variable kriter,charge
variable k,n, kc, nc, ratio, scale
nvar findindata
string lesprobas, lesmasses, ciblep,ciblem, tagg=""
variable/G lastincox, lastincoy
Killwaves/Z simu_proba, simu_mass
Make/O/D/N=1 simu_proba, simu_mass, simu_probadm
simu_proba=1
simu_mass=0
simu_probadm=defmass(simu_mass)	
wave/Z simu_proba, simu_mass, simu_probadm, masseel, transmass, transproba, incoy, incox//, probael
lesprobas=Wavelist("*_probabilites",";","") //recupere la liste des combinaison telles que genere par choisylecrible
lesmasses=Wavelist("*_massexact",";","")
n=itemsinlist(lesprobas)//les compte
	for(k=0;k<n;k+=1) //pour chaque élément
	ciblep = stringfromlist(k,lesprobas) //on regarde les probas des combinaisons isotopiques de cet element
	ciblem = stringfromlist(k,lesmasses)
	Duplicate/O $ciblep, probael
	Duplicate/O $ciblem, masseel
	nc=numpnts(masseel) //on compte le nombre de combinaisons pour l'element	
		for(kc=0;kc<nc;kc+=1) //pour chaque combinaison
		ratio=probael[kc]/wavemax(probael) //on evalue la probabilité par rapport à la probamax
			if (ratio>=kriter) //qui doit etre superieur ou egal au critère pour etre pris en compte
				tagg="augmasse"+num2str(k)+num2str(kc)
				Duplicate/O simu_mass, transmass
				transmass=transmass+masseel[kc]
				Duplicate/O transmass, $tagg
				Killwaves transmass
				tagg="augproba"+num2str(k)+num2str(kc)
				Duplicate/O simu_proba, transproba
				transproba=transproba*probael[kc]
				Duplicate/O transproba, $tagg
				Killwaves transproba
			endif
		endfor
		tagg=wavelist("augmasse*",";","")
		Concatenate/O/KILL/NP tagg, simu_mass
		tagg=wavelist("augproba*",";","")
		Concatenate/O/KILL/NP tagg, simu_proba
		Killwaves probael, masseel
	endfor
simu_mass= simu_mass-charge*0.00054857990946
	if (charge!=0)
		simu_mass= simu_mass/abs(charge)
	endif
Duplicate/O simu_mass deltam deltai test1//autoscale et calcule les ecarts en unité de masse ATTENTION PARAM ARBITRAIRE 0.5 et 0.8
wave/Z delta, wave0stripped, wave1stripped, test1
	switch(findindata)
		case 0:
			scale=wavemax(wave1stripped)/wavemax(simu_proba)
			break
		case 1:
			getaxis/W=elaborateur/Q right
			findvalue/V=(wavemax(simu_proba)) simu_proba
			//scale=wave1stripped[prochepic(simu_mass[V_value],0.1*V_max)]/wavemax(simu_proba)
			scale=incoy[0]/wavemax(simu_proba)
			break
		case 2:
			scale=incoy[0]/wavemax(simu_proba)
			break
	endswitch
simu_proba= simu_proba*scale
deltam=(simu_mass[x]-wave0stripped[prochepic(simu_mass[x],0.1*simu_proba[x])])
deltai=(simu_proba[x]-wave1stripped[prochepic(simu_mass[x],0.1*simu_proba[x])])
test1=sqrt(deltam[x]^2+deltai[x]^2)
labelize(simu_mass,simu_proba)
duplicate/O simu_mass, simu_probadm, simu_plot
simu_plot=1
simu_probadm=defmass(simu_mass)
lastincox=incox[0]
lastincoy=incoy[0]
end

Function SetKriter(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: 
		genestringmol()
		case 3: // Live update
		genestringmol()
			Variable dval = sva.dval/10
			String sval = sva.sval
			break
	endswitch

	return 0
End

Function Display_sel(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
NVAR dismode
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			labelize(simu_mass,simu_proba)
			break
	endswitch

	return 0
End

Window elaborateur() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(722.25,567.5,1071.75,704) wave1stripped vs wave0stripped
	AppendToGraph incoy vs incox
	AppendToGraph simu_proba vs simu_mass
	AppendToGraph simu_proba vs simu_mass
	AppendToGraph simu_profile_MixP4_outY vs simu_profile_MixP4_outX
	AppendToGraph simu_proba vs simu_mass
	ModifyGraph margin(top)=28
	ModifyGraph mode(wave1stripped)=8,mode(incoy)=3,mode(simu_proba)=3,mode(simu_proba#1)=3
	ModifyGraph mode(simu_proba#2)=8
	ModifyGraph marker(wave1stripped)=19,marker(incoy)=42,marker(simu_proba#2)=18
	ModifyGraph lSize(simu_profile_MixP4_outY)=2
	ModifyGraph rgb(wave1stripped)=(34816,34816,34816),rgb(incoy)=(36864,14592,58880)
	ModifyGraph rgb(simu_proba)=(65280,0,0),rgb(simu_proba#1)=(0,0,0),rgb(simu_profile_MixP4_outY)=(0,52224,52224)
	ModifyGraph rgb(simu_proba#2)=(0,52224,52224)
	ModifyGraph msize(wave1stripped)=1,msize(simu_proba)=3
	ModifyGraph mrkThick(wave1stripped)=1.1
	ModifyGraph opaque(wave1stripped)=1,opaque(incoy)=1
	ModifyGraph useMrkStrokeRGB(wave1stripped)=1,useMrkStrokeRGB(simu_proba#2)=1
	ModifyGraph textMarker(simu_proba)={labelmass,"default",3,0,1,0.00,15.00},textMarker(simu_proba#1)={labelpro,"default",1,0,1,0.00,2.00}
	ModifyGraph grid(left)=1
	ModifyGraph log(left)=1
	ModifyGraph nticks(bottom)=100
	ModifyGraph fStyle(left)=1
	ModifyGraph lblMargin(left)=116,lblMargin(bottom)=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2.42857,axOffset(bottom)=-1.33333
	ModifyGraph gridRGB(left)=(17408,17408,17408)
	ModifyGraph lblLatPos(left)=-1,lblLatPos(bottom)=369
	ModifyGraph axisEnab(bottom)={0.1,0.9}
	Label bottom "m/z"
	SetAxis left 1.65255451393238e-05,0.00165255451393238
	SetAxis bottom 409.86039875067,411.685298444511
	ErrorBars simu_proba#2 X,wave=(deltam,deltam)
	Legend/C/N=text0/J/D={1,0.5,0}/B=3/A=MC/X=48.35/Y=54.84 "\\s(simu_proba#2)Simu.\r\\s(simu_profile_MixP4_outY) Profile\r\\s(wave1stripped)Data"
	ControlBar 20
	SetVariable setvar0,pos={1,1},size={84,16},proc=SetVarProc_2,title="Width"
	SetVariable setvar0,limits={0,inf,0.1},value= simu_reso
	SetVariable setvar1,pos={86,1},size={93,16},proc=SetVarProc_2,title="Shape"
	SetVariable setvar1,limits={0,inf,1},value= simu_shape
	SetVariable setvar2,pos={180,1},size={84,16},proc=SetVarProc_2,title="Asym."
	SetVariable setvar2,limits={0,inf,1},value= simu_asym
	ValDisplay valdisp0,pos={264,1},size={99,14},title="FWHM"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp0,value= #"simu_prof_sigma*2.35482"
	PopupMenu popup0,pos={366,-1},size={73,21},proc=switch_profile
	PopupMenu popup0,mode=1,popvalue="Gaussian",value= #"\"Gaussian;Lorentz\""
EndMacro

function videagreg()
string lezmolz=wavelist("Formule_*",";",""),lamol
variable k,n=itemsinlist(lezmolz)
wave/T list_formules
for (k=0;k<n;k+=1)
	lamol=stringfromlist(k,lezmolz)
	removefromgraph/W=agregateur/Z $lamol
	Killwaves/Z $lamol
endfor
lezmolz=wavelist("DefmassForm_*",";","")
n=itemsinlist(lezmolz)
	for(k=0;k<n;k+=1)
		lamol=stringfromlist(k,lezmolz)
		removefromgraph/W=dmvm/Z $lamol
		Killwaves/Z $lamol
	endfor
lezmolz=wavelist("Molecule_*",";","")
n=itemsinlist(lezmolz)
	for(k=0;k<n;k+=1)
		lamol=stringfromlist(k,lezmolz)
		Killwaves/Z $lamol
	endfor
	for(k=n-1;k>-1;k-=1)
		MolbagWithdrawFrom({k},"Current")
	endfor
lezmolz=wavelist("Formule_*",";","")
n=itemsinlist(lezmolz)
	for(k=0;k<n;k+=1)
		lamol=stringfromlist(k,lezmolz)
		Killwaves/Z $lamol
	endfor	
Make/O/N=0/T List_formules
Make/D/O/N=0 Index_formules
make/O/D/N=(0,5) List_targets
//Setaxis/W=agregateur bottom wavemin(wave0stripped), wavemax(wave0stripped)
majliste()
//fullsavebag("Current")
Majlegende()
end

Window agregateur() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(331.5,42.5,819,484.25) wave1stripped vs wave0stripped
	AppendToGraph topbrighty vs topbrightx
	AppendToGraph treeintens vs tree
	AppendToGraph incoy vs incox
	AppendToGraph roi1 vs roi0
	AppendToGraph basey vs basex
	AppendToGraph simu_profile_MixP4_outY vs simu_profile_MixP4_outX
	AppendToGraph simu_proba vs simu_mass
	ModifyGraph gbRGB=(47872,47872,47872)
	ModifyGraph mode(wave1stripped)=8,mode(topbrighty)=3,mode(treeintens)=8,mode(incoy)=4
	ModifyGraph mode(roi1)=8,mode(basey)=3,mode(simu_proba)=8
	ModifyGraph marker(wave1stripped)=19,marker(topbrighty)=18,marker(treeintens)=55
	ModifyGraph marker(incoy)=42,marker(roi1)=19,marker(basey)=43,marker(simu_proba)=18
	ModifyGraph lSize(roi1)=2
	ModifyGraph rgb(wave1stripped)=(47872,47872,47872),rgb(topbrighty)=(0,0,0),rgb(incoy)=(36864,14592,58880)
	ModifyGraph rgb(roi1)=(36864,14592,58880),rgb(simu_profile_MixP4_outY)=(0,65280,65280)
	ModifyGraph rgb(simu_proba)=(0,52224,52224)
	ModifyGraph msize(wave1stripped)=1,msize(incoy)=5,msize(roi1)=3,msize(basey)=5
	ModifyGraph mrkThick(wave1stripped)=1,mrkThick(roi1)=1,mrkThick(simu_proba)=1
	ModifyGraph opaque(incoy)=1
	ModifyGraph hideTrace(basey)=1,hideTrace(simu_profile_MixP4_outY)=1
	ModifyGraph useMrkStrokeRGB(wave1stripped)=1,useMrkStrokeRGB(treeintens)=1,useMrkStrokeRGB(roi1)=1
	ModifyGraph useMrkStrokeRGB(simu_proba)=1
	ModifyGraph mrkStrokeRGB(simu_proba)=(65535,65535,65535)
	ModifyGraph zmrkSize(treeintens)={treeerror,0,10,6,0}
	ModifyGraph zColor(wave1stripped)={logwave1,-0.0772956579325359,*,YellowHot},zColor(treeintens)={treeColor,*,*,Web216,1}
	ModifyGraph zColorMin(wave1stripped)=(65535,65535,65535)
	ModifyGraph textMarker(topbrighty)={topbrightx,"default",0,90,1,0.00,0.00}
	ModifyGraph grid(left)=1
	ModifyGraph log(left)=1
	ModifyGraph zero(left)=1
	ModifyGraph nticks(bottom)=50
	ModifyGraph fStyle=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-1.85714,axOffset(bottom)=-1
	ModifyGraph axThick(left)=2
	ModifyGraph gridRGB(left)=(17408,17408,17408)
	ModifyGraph gridStyle(left)=1
	ModifyGraph gridHair(left)=3
	ModifyGraph lblPos(left)=45,lblPos(bottom)=44
	ModifyGraph tkLblRot(bottom)=90
	ModifyGraph axisEnab(bottom)={0.05,0.95}
	Label left "Intensity"
	Label bottom "m/z"
	SetAxis/A=2 left
	SetAxis bottom 150.002339,699.521952
	ControlBar 40
	Button button0,pos={1,2},size={80,20},proc=agregaReset,title="Zoom Reset"
	CheckBox check0,pos={1,24},size={61,14},proc=CheckProc,title="Logscale"
	CheckBox check0,variable= logscaleON
	SetVariable setvar0,pos={86,2},size={146,16},proc=SetVarProc,title="Top Bright Masses"
	SetVariable setvar0,limits={0,inf,1},value= top2dis
	ValDisplay valdisp0,pos={310,2},size={74,14},title="Y:",frame=0
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"pozint[0]"
	ValDisplay valdisp1,pos={389,2},size={72,14},title="X:",frame=0
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000},value= #"pozmasse[0]"
	ValDisplay valdisp2,pos={235,2},size={67,14},title="of"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp2,value= #"numpnts(wave1stripped)"
	PopupMenu popup0,pos={87,18},size={120,21},proc=switchpeakstick,title="Active mode"
	PopupMenu popup0,mode=1,popvalue="Sticks",value= #"\"Sticks;Profile\""
	SetWindow kwTopWin,hook=getagreg,hookevents=7
	SetWindow kwTopWin,hook(agreg)=getagreg
EndMacro

function majliste()
wave/T list_formules, label_formules
wave/Z index_formules, list_targets
variable n=numpnts(list_formules)
Make/O/T/N=4 listboxformulestitles={"","Name","Mass","dppm"}
make/O/T/N=(n,4) listboxformules=""
make/O/D/N=(n,4) selboxformules=0
listboxformules[][1]=list_formules[x]
listboxformules[][2]=num2str(index_formules[x])
listboxformules[][3]=num2str(list_targets[x][3])
selboxformules[][0]=32
insertpoints 0,1, listboxformules, selboxformules
selboxformules[0][]=32
//fullsavebag("Current")
end

function majlegende()
string lesmoles=wavelist("Formule_*",";",""),lamole
string lesdefmoles=wavelist("DefmassForm_*",";",""),ladefmole
variable k, n=itemsinlist(lesmoles)/2,i
wave/T listboxformules
wave/Z selboxformules
ValDisplay valdisp0 value=calcreqres(), win=molmanager
//Legend/W=agregateur/C/N=text0/J/F=0/A=MC ""
	for(k=0;k<n;k+=1)
		lamole=stringfromlist(2*k,lesmoles)
		ladefmole=stringfromlist(k,lesdefmoles)
		string labelstr="\JC"
			if(selboxformules[k+1][0]==48)
				for(i=1;i<4;i+=1)
					if(selboxformules[0][i]==48)
						labelstr+=listboxformules[k+1][i]+"\r"
					endif
				endfor
				labelstr=labelstr[0,strlen(labelstr)-2]
				findvalue/V=(wavemax($lamole)) $lamole
				Tag/W=agregateur/C/N=$lamole/TL=0/F=2/B=3/D={1,0.5,0} $lamole, V_value, labelstr
				Tag/W=dmvm/C/N=$ladefmole/TL=0/F=2/B=3/D={1,0.5,0} $ladefmole, V_value, labelstr
			else
				tag/W=agregateur/K/N=$lamole
				tag/W=dmvm/K/N=$ladefmole
			endif
	endfor
end


Function Vider(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			videagreg()
			break
	endswitch

	return 0
End

function seuilbruitFAT(rawdatay)
wave/Z rawdatay
Make/O/D/N=0 histo
wave/Z histo
Duplicate/O rawdatay logwave1
logwave1=log(rawdatay+1e-32)
histogram/B={wavemin(logwave1),0.1,10*(wavemax(logwave1)-wavemin(logwave1))} logwave1 histo
Killwaves/Z logwave1
differentiate histo
findpeak/Q/B=2/M=(wavemin(histo)/4)/N/R=[numpnts(histo)-1,0] histo
//Killwaves/Z histo
return V_PeakLoc
end

function deletenoiseFAT(rawy,rawx)
wave/Z rawy,rawx
duplicate/O rawy transy
duplicate/O rawx transx
sort rawy,transy,transx
variable seuil= 10^(seuilbruitFAT(rawy))
findlevel/Q transy seuil
deletepoints 0, v_levelx, transy,transx
sort transx,transy,transx
duplicate/O transy roi1
duplicate/O transx roi0
killwaves transy, transx
end

function keepeaks()
wave/Z wave1stripped,wave0stripped, w_findlevels, transy, transx
duplicate/O wave1stripped, wave1stripped_DIF
differentiate wave1stripped/X=wave0stripped/D=wave1stripped_DIF
findlevels/Q wave1stripped_DIF 0
w_findlevels/=(wave1stripped_DIF(floor(w_findlevels))>0)
wavetransform zapINFs w_findlevels
duplicate/O w_findlevels, transy, transx
transy=wave1stripped[w_findlevels]
transx=wave0stripped[w_findlevels]
duplicate/O transy, roi1
duplicate/O transx, roi0
killwaves transy, transx, wave1stripped_DIF
end

function seuilbruitFINE(rawdatay)
wave/Z rawdatay
Make/O/D/N=0 histo
wave/Z histo
Duplicate/O rawdatay logwave1
logwave1=log(rawdatay)
histogram/B={wavemin(logwave1),0.001,abs(wavemax(logwave1)-wavemin(logwave1))/0.001} logwave1 histo
findlevels/Q histo wavemax(histo)/2
Killwaves/Z logwave1, histo
return wavemax(w_findlevels)
end

function deletenoiseFINE(rawy,rawx)
wave/Z rawy,rawx
duplicate/O rawy transy
duplicate/O rawx transx
sort rawy,transy,transx
variable seuil= 10^(seuilbruitFINE(rawy))
findlevel/Q transy seuil
deletepoints 0, v_levelx, transy,transx
sort transx,transy,transx
duplicate/O transy roi1
duplicate/O transx roi0
killwaves transx,transy
end

function renomeactive(newname)
string newname
wave/Z wave1stripped,wave0stripped
string newnamey, newnamex
newnamey=newname+"_1stripped"
newnamex=newname+"_0stripped"
duplicate/O wave1stripped, $newnamey
duplicate/O wave0stripped, $newnamex
listelesdonnees()
end

Window panel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(7,797,942,978)
	Button H,pos={1,1},size={25,25},proc=elemH,title="H",fColor=(60928,60928,60928)
	Button He,pos={807,1},size={25,25},proc=elemHe,title="He"
	Button Li,pos={1,27},size={25,25},proc=elemLi,title="Li"
	Button Be,pos={27,27},size={25,25},proc=elemBe,title="Be"
	Button B,pos={677,27},size={25,25},proc=elemB,title="B"
	Button C,pos={703,28},size={25,25},proc=elemC,title="C",fColor=(4352,4352,4352)
	Button N,pos={729,27},size={25,25},proc=elemN,title="N"
	Button N,fColor=(48896,65280,57344)
	Button O,pos={755,27},size={25,25},proc=elemO,title="O"
	Button O,fColor=(65280,48896,48896)
	Button F,pos={781,27},size={25,25},proc=elemF,title="F"
	Button Ne,pos={807,27},size={25,25},proc=elemNe,title="Ne"
	Button Na,pos={1,53},size={25,25},proc=elemNa,title="Na"
	Button Mg,pos={27,53},size={25,25},proc=elemMg,title="Mg"
	Button Al,pos={677,53},size={25,25},proc=elemAl,title="Al"
	Button Si,pos={703,53},size={25,25},proc=elemSi,title="Si"
	Button P,pos={729,53},size={25,25},proc=elemP,title="P"
	Button S,pos={755,53},size={25,25},proc=elemS,title="S"
	Button S,fColor=(65280,65280,48896)
	Button Cl,pos={781,53},size={25,25},proc=elemCl,title="Cl"
	Button Ar,pos={807,53},size={25,25},proc=elemAr,title="Ar"
	Button K,pos={1,79},size={25,25},proc=elemK,title="K"
	Button Ca,pos={27,79},size={25,25},proc=elemCa,title="Ca"
	Button Sc,pos={417,79},size={25,25},proc=elemSc,title="Sc"
	Button Ti,pos={443,79},size={25,25},proc=elemTi,title="Ti"
	Button V,pos={469,79},size={25,25},proc=elemV,title="V"
	Button Cr,pos={495,79},size={25,25},proc=elemCr,title="Cr"
	Button Mn,pos={521,79},size={25,25},proc=elemMn,title="Mn"
	Button Fe,pos={547,79},size={25,25},proc=elemFe,title="Fe"
	Button Co,pos={573,79},size={25,25},proc=elemCo,title="Co"
	Button Ni,pos={599,79},size={25,25},proc=elemNi,title="Ni"
	Button Cu,pos={625,79},size={25,25},proc=elemCu,title="Cu"
	Button Zn,pos={651,79},size={25,25},proc=elemZn,title="Zn"
	Button Ga,pos={677,79},size={25,25},proc=elemGa,title="Ga"
	Button Ge,pos={703,79},size={25,25},proc=elemGe,title="Ge"
	Button As,pos={729,79},size={25,25},proc=elemAs,title="As"
	Button Se,pos={755,79},size={25,25},proc=elemSe,title="Se"
	Button Br,pos={781,79},size={25,25},proc=elemBr,title="Br"
	Button Kr,pos={807,79},size={25,25},proc=elemKr,title="Kr"
	Button Rb,pos={1,105},size={25,25},proc=elemRb,title="Rb"
	Button Sr,pos={27,105},size={25,25},proc=elemSr,title="Sr"
	Button Y,pos={417,105},size={25,25},proc=elemY,title="Y"
	Button Zr,pos={443,105},size={25,25},proc=elemZr,title="Zr"
	Button Nb,pos={469,105},size={25,25},proc=elemNb,title="Nb"
	Button Mo,pos={495,105},size={25,25},proc=elemMo,title="Mo"
	Button Tc,pos={521,105},size={25,25},proc=elemTc,title="Tc"
	Button Ru,pos={547,105},size={25,25},proc=elemRu,title="Ru"
	Button Rh,pos={573,105},size={25,25},proc=elemRh,title="Rh"
	Button Pd,pos={599,105},size={25,25},proc=elemPd,title="Pd"
	Button Ag,pos={625,105},size={25,25},proc=elemAg,title="Ag"
	Button Cd,pos={651,105},size={25,25},proc=elemCd,title="Cd"
	Button In,pos={677,105},size={25,25},proc=elemIn,title="In"
	Button Sn,pos={703,105},size={25,25},proc=elemSn,title="Sn"
	Button Sb,pos={729,105},size={25,25},proc=elemSb,title="Sb"
	Button Te,pos={755,105},size={25,25},proc=elemTe,title="Te"
	Button I,pos={781,105},size={25,25},proc=elemI,title="I"
	Button Xe,pos={807,105},size={25,25},proc=elemXe,title="Xe"
	Button Cs,pos={1,131},size={25,25},proc=elemCs,title="Cs"
	Button Ba,pos={27,131},size={25,25},proc=elemBa,title="Ba"
	Button La,pos={53,131},size={25,25},proc=elemLa,title="La"
	Button Ce,pos={79,131},size={25,25},proc=elemCe,title="Ce"
	Button Pr,pos={105,131},size={25,25},proc=elemPr,title="Pr"
	Button Nd,pos={131,131},size={25,25},proc=elemNd,title="Nd"
	Button Pm,pos={157,131},size={25,25},proc=elemPm,title="Pm"
	Button Sm,pos={183,131},size={25,25},proc=elemSm,title="Sm"
	Button Eu,pos={209,131},size={25,25},proc=elemEu,title="Eu"
	Button Gd,pos={235,131},size={25,25},proc=elemGd,title="Gd"
	Button Tb,pos={261,131},size={25,25},proc=elemTb,title="Tb"
	Button Dy,pos={287,131},size={25,25},proc=elemDy,title="Dy"
	Button Ho,pos={313,131},size={25,25},proc=elemHo,title="Ho"
	Button Er,pos={339,131},size={25,25},proc=elemEr,title="Er"
	Button Tm,pos={365,131},size={25,25},proc=elemTm,title="Tm"
	Button Yb,pos={391,131},size={25,25},proc=elemYb,title="Yb"
	Button Lu,pos={417,131},size={25,25},proc=elemLu,title="Lu"
	Button Hf,pos={443,131},size={25,25},proc=elemHf,title="Hf"
	Button Ta,pos={469,131},size={25,25},proc=elemTa,title="Ta"
	Button W,pos={495,131},size={25,25},proc=elemW,title="W"
	Button Re,pos={521,131},size={25,25},proc=elemRe,title="Re"
	Button Os,pos={547,131},size={25,25},proc=elemOs,title="Os"
	Button Ir,pos={573,131},size={25,25},proc=elemIr,title="Ir"
	Button Pt,pos={599,131},size={25,25},proc=elemPt,title="Pt"
	Button Au,pos={625,131},size={25,25},proc=elemAu,title="Au"
	Button Hg,pos={651,131},size={25,25},proc=elemHg,title="Hg"
	Button Tl,pos={677,131},size={25,25},proc=elemTl,title="Tl"
	Button Pb,pos={703,131},size={25,25},proc=elemPb,title="Pb"
	Button Bi,pos={729,131},size={25,25},proc=elemBi,title="Bi"
	Button Po,pos={755,131},size={25,25},proc=elemPo,title="Po"
	Button At,pos={781,131},size={25,25},proc=elemAt,title="At"
	Button Rn,pos={807,131},size={25,25},proc=elemRn,title="Rn"
	Button Fr,pos={1,157},size={25,25},proc=elemFr,title="Fr"
	Button Ra,pos={27,157},size={25,25},proc=elemRa,title="Ra"
	Button Ac,pos={53,157},size={25,25},proc=elemAc,title="Ac"
	Button Th,pos={79,157},size={25,25},proc=elemTh,title="Th"
	Button Pa,pos={105,157},size={25,25},proc=elemPa,title="Pa"
	Button U,pos={131,157},size={25,25},proc=elemU,title="U"
	Button Np,pos={157,157},size={25,25},proc=elemNp,title="Np"
	Button Pu,pos={183,157},size={25,25},proc=elemPu,title="Pu"
	Button Am,pos={209,157},size={25,25},proc=elemAm,title="Am"
	Button Cm,pos={235,157},size={25,25},proc=elemCm,title="Cm"
	Button Bk,pos={261,157},size={25,25},proc=elemBk,title="Bk"
	Button Cf,pos={287,157},size={25,25},proc=elemCf,title="Cf"
	Button Es,pos={313,157},size={25,25},proc=elemEs,title="Es"
	Button Fm,pos={339,157},size={25,25},proc=elemFm,title="Fm"
	Button Md,pos={365,157},size={25,25},proc=elemMd,title="Md"
	Button No,pos={391,157},size={25,25},proc=elemNo,title="No"
	Button Lr,pos={417,157},size={25,25},proc=elemLr,title="Lr"
	Button Rf,pos={443,157},size={25,25},proc=elemRf,title="Rf"
	Button Db,pos={469,157},size={25,25},proc=elemDb,title="Db"
	Button Sg,pos={495,157},size={25,25},proc=elemSg,title="Sg"
	Button Bh,pos={521,157},size={25,25},proc=elemBh,title="Bh"
	Button Hs,pos={547,157},size={25,25},proc=elemHs,title="Hs"
	Button Mt,pos={573,157},size={25,25},proc=elemMt,title="Mt"
	Button Uun,pos={599,157},size={25,25},proc=elemUun,title="Uun"
	Button Uuu,pos={625,157},size={25,25},proc=elemUuu,title="Uuu"
	Button Uub,pos={651,157},size={25,25},proc=elemUub,title="Uub"
	Button Uuq,pos={677,157},size={25,25},proc=elemUuq,title="Uuq"
	Button Uuh,pos={703,157},size={25,25},proc=elemUuh,title="Uuh"
	Button Zeromol,pos={833,53},size={102,50},proc=ButtonReset,title="Reset"
	Button Zeromol,help={"Clears current Molecule"},fColor=(65280,43520,0)
	SetVariable setvarkrikri,pos={84,53},size={269,16},bodyWidth=45,proc=SetKriter,title="Isotopic combo probability to take into account"
	SetVariable setvarkrikri,help={"Sets threshold of isotopical combination probability of each element to be reached to be taken into account in the simulation"}
	SetVariable setvarkrikri,limits={0,1,0.01},value= krikriter
	CheckBox dis_selector,pos={53,79},size={147,14},proc=Display_sel,title="Display relative to maximum"
	CheckBox dis_selector,help={"Test me"},variable= dismode
	Button adder,pos={833,1},size={102,51},proc=Ajouter,title="Add to agregator"
	Button adder,help={"Sends the current Molecule to the agregator window that stores selected isotopical simulations"}
	Button adder,fColor=(32768,65280,0)
	SetVariable charge_setter,pos={60,105},size={91,16},bodyWidth=36,proc=RegleCharge,title="Ion charge"
	SetVariable charge_setter,help={"Sets the current ion charge to work with. Impacts calculation and mass simulation."}
	SetVariable charge_setter,limits={-5,5,1},value= charge
	Button toperm,pos={53,1},size={102,51},proc=ButtonAddtoPerm,title="Add as factor"
	Button toperm,help={"Sends the current Molecule (right) to the factors list"}
	Button toperm,fColor=(32768,65280,0)
	TitleBox title0,pos={313,1},size={21,25}
	TitleBox title0,help={"Displays the current Molecule. It can come from anywhere and go to anywhere. Simply click on elements button to modify. Isotopical combo is simulated and displayed above."}
	TitleBox title0,labelBack=(65535,65535,65535),variable= nommol
	Button Recalibrator,pos={833,105},size={51,51},proc=ButtonCalib,title="Recal"
	Button Recalibrator,help={"Shifts the current sample data in mass so the closest peak to the current formula is set at the calculated mass."}
	Button Recalibrator,fColor=(65280,0,0)
	Button Recalibrator1,pos={884,105},size={51,51},proc=ButtonHardCalib,title="Hard\rRecal"
	Button Recalibrator1,help={"Shifts the current sample data in mass so the closest peak to the current formula is set at the calculated mass."}
	Button Recalibrator1,fColor=(65280,0,0)
	Button msmser,pos={833,157},size={102,25},proc=gotoMSMS,title="doMSMS"
	Button msmser,fColor=(57344,65280,48896)
EndMacro


function effacemolecule()
string formuletokill, formulemasstokill, moleculetokill, defmasstokill
wave/Z list_formules, index_formules//, List_targets
controlinfo /W=molmanager list0
variable n
	if (V_value>0)
		formuletokill=stringfromlist(2*(V_value-1),wavelist("Formule_*",";",""))
		formulemasstokill=stringfromlist(2*(V_value-1)+1,wavelist("Formule_*",";",""))
		moleculetokill=stringfromlist(V_value-1,wavelist("Molecule_*",";",""))
		n=itemsinlist(wavelist("Molecule_*",";",""))
		defmasstokill=stringfromlist(V_value-1,wavelist("DefmassForm_*",";",""))
		Removefromgraph/W=agregateur $formuletokill
		Removefromgraph/W=dmvm $defmasstokill
		Killwaves/Z $formuletokill,$formulemasstokill, $moleculetokill, $defmasstokill
		deletepoints/M=0 V_value-1, 1, list_formules, index_formules, List_targets
		MolbagWithdrawFrom({(v_value-1)},"current")
		if(numpnts(list_targets)==0)
			make/O/D/N=(0,5) List_targets
		endif
	endif
	majformulewavesnames()
	majliste()
	majlegende()
end

function majformulewavesnames()
string lesformules=wavelist("Formule_*",";",""), formuletorename, nouveaublaz, formuletorenamemass, nouveaublazmass, moleculetorename, newmolname, defmasstorename, newdefmass
string lesmol=wavelist("molecule_*",";",""), lesdef=wavelist("DefmassForm_*",";","")
variable k, n=itemsinlist(lesformules)
	for(k=0;k<n;k+=2)
		formuletorename=stringfromlist(k,lesformules)
		nouveaublaz="Formule_"+num2str(k/2+1)
		formuletorenamemass=stringfromlist(k+1,lesformules)
		nouveaublazmass="Formule_"+num2str(k/2+1)+"mass"
		moleculetorename=stringfromlist(k/2,lesmol)
		newmolname="Molecule_"+num2str(k/2+1)
		defmasstorename=stringfromlist(k/2,lesdef)
		newdefmass="DefmassForm_"+num2str(k/2+1)
			if (cmpstr(formuletorename,nouveaublaz)!=0)
				rename $formuletorename $nouveaublaz
				rename $formuletorenamemass $nouveaublazmass
				rename $moleculetorename $newmolname
				rename $defmasstorename $newdefmass
			endif
	endfor
end

Function DeleteMol(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave/Z list_formules
	switch( ba.eventCode )
		case 2: // mouse up
		controlinfo /W=molmanager list0
		if(numpnts(list_formules)!=0 && v_value>0)
			effacemolecule()
		endif
			break
	endswitch

	return 0
End

Function RegleCharge(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		genestringmol()
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			genestringmol()
			break
	endswitch

	return 0
End

function focusmol()
string masstofocus
controlinfo /W=molmanager list0
	if (V_value>0 && dimsize($S_value,0)>1)
		masstofocus=stringfromlist(2*(V_value-1)+1,wavelist("Formule_*",";",""))
		Setaxis/W=agregateur bottom wavemin($masstofocus)-2, wavemax($masstofocus)+2
		Setaxis/W=elaborateur bottom wavemin($masstofocus)-2, wavemax($masstofocus)+2
	else
	Setaxis/W=agregateur bottom wavemin(wave0stripped), wavemax(wave0stripped)
	endif

end


Window panel0() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1113,493,1267,979)
	ListBox list0,pos={2,1},size={153,334},proc=ListBoxProcAgregToElab
	ListBox list0,help={"Displays the list of formulas shown in the agregator. Click to select one."}
	ListBox list0,labelBack=(0,65280,0),fSize=12,frame=4,fStyle=0
	ListBox list0,listWave=root:List_formules,mode= 1,selRow= 0,special= {0,18,0}
	Button Delete,pos={1,436},size={76,51},proc=DeleteMol,title="Erase"
	Button Delete,help={"Removes selected formula from list and from agregator window"}
	Button Delete,fColor=(65280,0,0)
	Button centre,pos={2,335},size={76,51},proc=centrevue,title="Focus"
	Button centre,help={"Changes the agregator window settings to show the entire isotopical simulation. Focus when none selected restore default settings."}
	Button centre,fColor=(32768,65280,0)
	Button centre1,pos={77,335},size={76,51},proc=Switchmass,title="Check\r Isotopes\r Pattern"
	Button centre1,help={"Spans each peak of an isotopical simulation in the agregator window."}
	Button centre1,fColor=(32768,65280,0)
	Button autoscale,pos={1,385},size={76,51},proc=rescaleselected,title="Rescale\r selected"
	Button autoscale,help={"Sets the maximum of selected isotopical simulation in the agregator equal to the height of the closest peak in the dataset"}
	Button autoscale,fColor=(65280,43520,0)
	Button rescaleALL,pos={77,385},size={76,51},proc=rescaleall,title="Rescale all"
	Button rescaleALL,help={"Sets the maximum of each isotopical simulation in the agregator equal to the height of the closest peak in the dataset"}
	Button rescaleALL,fColor=(65280,43520,0)
	Button Emptier6,pos={77,436},size={76,51},proc=Vider,title="Empty\r agregator"
	Button Emptier6,help={"Removes every isotopical simulation from the agregator (top above)"}
	Button Emptier6,fColor=(65280,0,0)
EndMacro

Function centrevue(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			focusmol()
			break
	endswitch

	return 0
End

function prochepic(pos,seuil)//pos en masse mais réponse en point
variable pos,seuil
wave/Z wave1stripped, wave0stripped
findlevel/Q wave0stripped pos
return round(v_levelx)
end

function procheETintense(pos)
variable pos
wave/Z wave0stripped, dataindex
nvar charge
pos=round(charge*pos)/charge
duplicate/O dataindex transclass
transclass=abs(round(charge*wave0stripped(dataindex))/charge - pos)
findlevel/Q transclass 0
killwaves/Z transclass
return wave0stripped(dataindex[v_levelx])
end


function ANCIENprochepic(pos,seuil)//pos en masse mais réponse en point
variable pos,seuil
variable res=(wavemax(wave0stripped)-wavemin(wave0stripped))/numpnts(wave0stripped), agau, adro
wave/Z wave1stripped, wave0stripped
findvalue/V=(pos)/T=(2*res) wave0stripped
findpeak/B=(10*res)/M=(seuil)/Q/R=[v_value,0] wave1stripped
	if(V_flag==0)
		agau=V_peakloc
	else
		agau=0
	endif
findpeak/B=(10*res)/M=(seuil)/Q/R=[v_value,numpnts(wave1stripped)] wave1stripped
	if(V_flag==0)
		adro=v_peakloc
	else
		adro=numpnts(wave1stripped)
	endif
	
	if(abs(wave0stripped[agau]-pos)<abs(wave0stripped[adro]-pos))
		return agau
	else
		return adro
	endif
end

Function Switchmass(ba) : ButtonControl
	STRUCT WMButtonAction &ba
		nvar Cfocus

	switch( ba.eventCode )
		case 2: // mouse up
			Cfocus+=1
			isopattern()
			break
	endswitch

	return 0
End


function isopattern()
nvar Cfocus
string masstofocus
controlinfo /W=molmanager list0
	if (V_value>0 && dimsize($S_value,0)>1)
		masstofocus=stringfromlist(2*(V_value-1)+1,wavelist("Formule_*",";",""))
		Duplicate/O $masstofocus lamasse
		variable n=numpnts(lamasse)
			if(Cfocus>n-1)
				Cfocus=0
			endif
	endif
Setaxis/W=agregateur bottom lamasse[Cfocus]-0.05, lamasse[Cfocus]+0.05
end



Window panel1() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(960,797,1096,979)
	ListBox list0,pos={1,1},size={135,128},proc=choisyladata
	ListBox list0,help={"Displays available datasets. Click to change and refresh top ten peaks."}
	ListBox list0,listWave=root:liste_donnees,mode= 1,selRow= 1
	Button buttonLoadData,pos={1,130},size={35,25},proc=ButtonLoader,title="Load"
	Button buttonLoadData,help={"Load an XY text file with mass and intensities."}
	Button buttonLoadData,fColor=(32768,65280,0)
	Button buttonKeepeaks,pos={38,147},size={60,16},proc=ButtonPeaks,title="Peaks"
	Button buttonKeepeaks,help={"Find the centroid of each peak found in the current sample data by derivation."}
	Button buttonKeepeaks,fColor=(65280,43520,0)
	Button buttonKiller,pos={100,155},size={35,25},proc=ButtonKill,title="Kill"
	Button buttonKiller,help={"Removes the selected sample data."}
	Button buttonKiller,fColor=(65280,0,0)
	Button buttonFatNoise,pos={38,130},size={60,16},proc=ButtonKillFatNoise,title="Fat Noise"
	Button buttonFatNoise,help={"Remove the less intense points so the new intensity distribution is strictly convex"}
	Button buttonFatNoise,fColor=(65280,43520,0)
	Button buttonFineNoise,pos={38,164},size={60,16},proc=ButtonKillFineNoise,title="Fine Cut"
	Button buttonFineNoise,help={"Deletes the less intense points so the most intense peaks becomes twice less frequent."}
	Button buttonFineNoise,fColor=(65280,43520,0)
	Button buttonRename,pos={100,130},size={35,25},proc=ButtonDupeRname,title="Dupe"
	Button buttonRename,help={"Generates an other sample data from the current data."}
	Button buttonRename,fColor=(32768,65280,0)
	Button buttonNormalize,pos={1,155},size={35,25},proc=ButtonNorm,title="Norm"
	Button buttonNormalize,help={"Divides the current sample data by its maximum intensity"}
	Button buttonNormalize,fColor=(32768,65280,0)
EndMacro

threadsafe function defmass(mass)
variable mass
return mass-round(mass)
end

function intenstat(percent)
variable percent
wave/Z wave1stripped, dataindex
variable xeuma=sum(wave1stripped)
if(percent==1)
	return numpnts(wave1stripped)
else
	duplicate/O/FREE wave1stripped trans1
	trans1=wave1stripped[dataindex[x]]
	Integrate/DIM=0  trans1
	trans1/=xeuma
	findlevel/Q trans1 percent
	if(numtype(ceil(v_levelx)+1)!=0)
		return 1
	endif
return round(v_levelx)+1
endif
end

function rescaleformules()
variable k=0, n=numpnts(index_formules),maxi,scale
wave wave1stripped, wave0stripped, list_targets, index_formules, lawave
string formuletorescale
	for(k=0;k<n;k+=1)
		formuletorescale="Formule_"+num2str(k+1)
		maxi=wavemax($formuletorescale)
		duplicate/O $formuletorescale, lawave
		scale=wave1stripped[prochepic(index_formules[k],0.1*maxi)]
		list_targets[k][1]=wave0stripped[prochepic(index_formules[k],0.1*maxi)]
		list_targets[k][2]=wave1stripped[prochepic(index_formules[k],0.1*maxi)]
		list_targets[k][3]=1e6*(index_formules[k]-list_targets[k][1])/index_formules[k]
		lawave*=(scale/maxi)
		duplicate/O lawave, $formuletorescale
	endfor
killwaves lawave
majliste()
fullsavebag("Current")
majlegende()
end

function rescaleTHATformule()
controlinfo/W=molmanager list0
variable k=V_value-1, maxi,scale
wave wave1stripped, wave0stripped, list_targets, index_formules, lawave
string formuletorescale
if(k>-1)
formuletorescale="Formule_"+num2str(k+1)
maxi=wavemax($formuletorescale)
duplicate/O $formuletorescale, lawave
scale=wave1stripped[prochepic(index_formules[k],0.1*maxi)]
list_targets[k][1]=wave0stripped[prochepic(index_formules[k],0.1*maxi)]
list_targets[k][2]=wave1stripped[prochepic(index_formules[k],0.1*maxi)]
list_targets[k][3]=1e6*(index_formules[k]-list_targets[k][1])/index_formules[k]
lawave*=(scale/maxi)
duplicate/O lawave, $formuletorescale
killwaves lawave
endif
majliste()
fullsavebag("Current")
majlegende()
end

Function rescaleselected(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			Execute "rescaleTHATformule()"
			break
	endswitch

	return 0
End

Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			execute "rescaleformules()"
			break
	endswitch

	return 0
End

Function rescaleall(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			execute "rescaleformules()"
			break
	endswitch

	return 0
End

function addMOLtoperm()//ajoute une ligne à list_perm et crée une nouvelle molperm
svar nommol
wave/Z Molecule
wave/T list_molperm
variable n=itemsinlist(wavelist("molperm_*",";","")), len
string newname="Molperm_"+num2str(n)
Duplicate Molecule, $newname
insertpoints n, 1, list_molperm

len = strsearch(nommol,"+",0)//eliminons la charge dans le nom
	if (len != -1)
		list_molperm[n]=nommol[0,len-2]
	else
		len = strsearch(nommol,"-",0)
		if (len != -1)
			list_molperm[n]=nommol[0,len-2]
		else
			list_molperm[n]=nommol
		endif
	endif
calcperm()
end

function Molbag2perm(str)
string str
string molbag="Molbag_"+str
string listbag="Listbag_"+str
string newname, lenom
wave transbag=$molbag
wave/T translist=$listbag, list_molperm
emptyperm()
duplicate/O/T translist list_molperm
duplicate/O molecule transmol
variable n=numpnts(translist), k, len
for(k=0;k<n;k+=1)
	newname="Molperm_"+num2str(k)
	transmol=transbag[x][k]
	duplicate transmol $newname
	lenom=translist[k]
	len = strsearch(lenom,"+",0)//eliminons la charge dans le nom
	if (len != -1)
		list_molperm[k]=lenom[0,len-2]
	else
		len = strsearch(lenom,"-",0)
		if (len != -1)
			list_molperm[k]=lenom[0,len-2]
		else
			list_molperm[k]=lenom
		endif
	endif
endfor
calcperm()
killwaves/Z transmol
end

function calcperm()//genere et calcule perm
wave/T list_molperm
wave/Z transmol, mendev_masses
variable n=numpnts(list_molperm), k
string lesmol=wavelist("molperm_*",";",""), lamol, newmol
Make/O/D/N=(numpnts(mendev_masses)) transmass
wave/Z transmass
transmass=mendev_masses[x][0]+0 //se cantone aux isotopes majoritaires
Make/O/D/N=(n) perm
	for(k=0;k<n;k+=1) //pour chaque molperm présente, calcule la masse de la stoechiometrie en isotope majoritaire.
		lamol="molperm_"+num2str(k)
		duplicate/O $lamol, transmol
		transmol*=transmass
		perm[k]=sum(transmol)
	endfor
Duplicate/O perm transindex
Makeindex perm transindex
rotate -1, transindex
duplicate/O perm transperm
Make/O/T/N=(n) translist
transperm=perm(transindex[x])
translist=list_molperm(transindex[x])
	for(k=0;k<n;k+=1)
		lamol="molperm_"+num2str(transindex[k])
		newmol="transmol_"+num2str(k)
		rename $lamol $newmol
	endfor
	for(k=0;k<n;k+=1)
		lamol="transmol_"+num2str(k)
		newmol="molperm_"+num2str(k)
		rename $lamol $newmol
	endfor
duplicate/O/T translist list_molperm
duplicate/O transperm perm
Make/O/T/N=(n,2) minmax
Make/O/D/N=(n,2) driveminmax
minmax[][0]="0"
minmax[][1]="100"
minmax[n-1][0]="0"
minmax[n-1][1]="1000"
driveminmax=2
driveminmax[n-1][]=0
if(n!=0)
concatenate/O sortlist(wavelist("Molperm_*",";","")), Molbag_Factors
else
Make/O/N=(114,0) Molbag_Factors
endif
duplicate/O list_molperm Listbag_Factors
Make/O/D/N=(n) Chargebag_Factors=0
Make/O/D/N=(n,5) Ciblebag_Factors=0
Ciblebag_Factors[][1]=perm[x]
Ciblebag_Factors[][2]=1
killwaves translist, transperm, transindex, transmol, transmass
end

function killthismol(this)
variable this
string lamol="molperm_"+num2str(this), lesmol=wavelist("molperm_*",";",""), newmol
variable, k, n=itemsinlist(lesmol)
killwaves $lamol
deletepoints this,1, perm, list_molperm
	for(k=this+1;k<n;k+=1)
		lamol="molperm_"+num2str(k)
		newmol="molperm_"+num2str(k-1)
		rename $lamol $newmol
	endfor
end

function permmax(inco)
variable inco
wave/Z perm
duplicate/O perm, maxperm
wave/Z maxperm
maxperm= ceil(inco/perm)
end

function trouvemass(inco)
variable inco
nvar charge
wave/Z perm, maxperm,masstosoustract//,propindex
variable k,n=numpnts(perm), j=0, i=1, sab=ticks, m, a
string/G compo=""
Make/O/D/N=(n) prod //Calcul le produit cumulatif du nombre max de chaque element
prod=1 // initialise
for(k=1;k<n;k+=1) //pour chaque masse au dela de la première
	prod[k]=prod[k-1]*(maxperm[k-1]+1) //calcul le produit des k-1 premières masses permises, le "stretch"
endfor
i=prod[n-1] //i est le nombre de ligne
killwaves/Z coefz, massetot //detruire les waves est plus rapide que de les annuler
Make/O/W/N=(0,n) coefz //genere une ligne de coef
Make/O/D/N=(0) massetot //genere la ligne de masse qui deviendra celle des defauts de masse
Make/O/D/N=(1,n-1) con //genere une ligne de transfert qui sert à plusieurs choses
for(j=0,a=0;j<i;j+=1,a+=1) //-----Début du calcul
	insertpoints a,1, coefz,massetot
	coefz[a][]=mod(trunc(j/prod[y]),maxperm[y]+1) //rempli la ligne de coef en garantissant le crible exhaustif sur les masses sauf la plus faible (dioph)
	con=coefz[a][y]*perm[y] //calcul la masse ajouté par élément
	massetot[a]=sum(con) //calcul de la masse
	coefz[a][n-1]=max(round((inco-massetot[a])/perm[n-1]),0) //libère un degré de liberté (la masse la plus faible sert à combler)
	massetot[a]+=coefz[a][n-1]*perm[n-1]//ajuste la masse
	if (abs(inco-massetot[a])>=perm[n-1]/2) //ATTENTION 0.5 PARAMETRE ARBITRAIRE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//si la masse calculée depasse on prend des mesures
		for(m=n-2;inco>sum(con,m,n-2) && m>0;m=m-1)	 //on cherche à partir de quel élément le dépassement de la masse est garanti
		endfor
		j=prod[m+1]*ceil(j/prod[m+1]) //on saute ici pour aller au j où le nombre d'élément dont le rajout garantissait le depassement retombe à zéro
		coefz[a][]=mod(trunc(j/prod[y]),maxperm[y]+1) // on refait l'opération de calcul pour ne pas perdre de temps et ne pas retomber surune ligne tocarde
		con=coefz[a][y]*perm[y]
		massetot[a]=sum(con)
		coefz[a][n-1]=max(round((inco-massetot[a])/perm[n-1]),0)
		massetot[a]+=coefz[a][n-1]*perm[n-1]
	endif
endfor
duplicate/o massetot lesmassescriblees
variable offset=- charge*0.00054857990946 +sum(masstosoustract)
variable scale=max(1,abs(charge))
fastop lesmassescriblees=(offset)+(scale)*lesmassescriblees//lesmassescriblees*max(1,abs(charge)) - charge*0.00054857990946 +sum(masstosoustract)
multithread massetot=(massetot-inco)*1e6/(massetot+sum(masstosoustract)) //on transforme en défaut de masse
duplicate/O massetot propindex 
multithread propindex=abs(propindex)
makeindex propindex propindex
return (ticks-sab)*(1/60)
end

function generelistpropOBSO_SINCE_2013(nb)// nb est le nombre de propositions souhaitées
variable nb
svar nommol
nvar charge
variable killen=numpnts(propindex)-nb, i=numpnts(perm), k, j
string newmolpropname, transname
wave/Z coefz, propindex, molecule, massetot,lesmassescriblees
wave/T minmax
deletepoints nb, killen, propindex
Make/O/T/N=(nb,2) list_prop
list_prop=""
Duplicate/O Molecule, transmolprop, newmolprop
	for(k=0;k<nb;k+=1) //iteration sur les propositions
		newmolpropname="molprop_"+num2str(k)
		newmolprop=0
		transmolprop=0
		for(j=0;j<i;j+=1)
			transname="molperm_"+num2str(j)
			Duplicate/O $transname, transmolprop
			newmolprop+=(str2num(minmax[j][0])+coefz[propindex[k]][j])*transmolprop
		endfor
		Duplicate/O newmolprop, molecule, $newmolpropname
		genestringmol()
		list_prop[k][0]=nommol
		molecule=0
		tuelesspec()
	endfor
	genestringmol()
	labelize(simu_mass,simu_proba)
	duplicate/o propindex masseclassee
	masseclassee=lesmassescriblees[propindex[x]]
	propindex=massetot[propindex[x]]
	Duplicate/O propindex deltappm
	list_prop[][1]=num2str(deltappm[x])
	killwaves transmolprop, newmolprop
end

function generelistprop(nb)// nb est le nombre de propositions souhaitées
variable nb
svar nommol
nvar charge
variable killen=numpnts(propindex)-nb, i=numpnts(perm), k, j, t=ticks
string newmolpropname, transname
wave/Z coefz, propindex, molecule, massetot,incox,incoy
wave/T minmax
//deletepoints nb, killen, propindex                                     // roland a mis ceci en commentaire pour avoir de la réserve de proposition dans le patch qui suit
Make/O/T/N=(nb,2) list_prop
list_prop=""
//ajout 14 fevrier 2019 pour modernisation doubloneuse
Make/O/D/N=(i)/FREE minimum=str2num(minmax[x][0])
Make/O/D/N=(i)/FREE zfeoc=0
////////
Duplicate/O Molecule, transmolprop, newmolprop
                for(k=0;k<nb;k+=1) //iteration sur les propositions
                               newmolpropname="molprop_"+num2str(k)
                               newmolprop=0
                               transmolprop=0
                               for(j=0;j<i;j+=1)
                                               transname="molperm_"+num2str(j)
                                               Duplicate/O $transname, transmolprop
                                               newmolprop+=(minimum[j]+coefz[propindex[k]][j])*transmolprop
                               endfor
                               zfeoc=minimum[x]+coefz[propindex[k]][x]
                               Duplicate/O/FREE weird2canon(zfeoc,"Factors"), fullprop
                               Duplicate/O newmolprop, molecule, $newmolpropname
                               //genestringmol() //commenté le 4 fevrier2019 pour obsolescence
                               nommol=mol2str(molecule,charge)
                               //// patch roland pour eviter les propositions identiques fevrier 2013
                               if (k>0)
                               			if (cmpstr(list_prop[k-1][0], nommol)==0)
                                                               nb-=1
                                                               deletepoints/M=0 k, 1, propindex, list_prop
                                                                k-=1
                                               else
                                                               list_prop[k][0]=nommol
                                                               MolbagAddTo(molecule,charge,nommol,{0,incox[0],incoy[0],massetot[propindex[k]],0},"Propositions")
                                               endif
                               else
                                               list_prop[k][0]=nommol
                                               MolbagAddTo(molecule,charge,nommol,{0,incox[0],incoy[0],massetot[propindex[k]],0},"Propositions")
                               endif
                               /////fin patch roland
                               //molecule=0 //commenté le 4 fevrier2019 pour obsolescence
                               //tuelesspec() //commenté le 4 fevrier2019 pour obsolescence
                endfor
                //genestringmol() //commenté le 4 fevrier2019 pour obsolescence
                //labelize(simu_mass,simu_proba) //commenté le 4 fevrier2019 pour obsolescence
                propindex=massetot[propindex[x]]
                Duplicate/O propindex deltappm
                list_prop[][1]=num2str(deltappm[x])
                killwaves transmolprop, newmolprop, propindex
//print ticks-t
end

function razprop()
string lesmolprop=wavelist("molprop_*",";",""), molproptokill
wave/Z deltappm
wave/T list_prop
variable k, n=itemsinlist(lesmolprop)
killwaves/Z deltappm
	for(k=0;k<n;k+=1)
		molproptokill=stringfromlist(k,lesmolprop)
		killwaves $molproptokill
	endfor
	Make/O/T/N=(0,2) list_prop
if(waveexists(Ciblebag_Propositions))
	n=dimsize(Ciblebag_Propositions,0)
	make/FREE/O/D/N=(n) mol2kill=x
	MolbagWithdrawFrom(mol2kill,"Propositions")
endif
end

macro nantozero(onde)
string onde
variable k=0
variable n=numpnts($onde)
do
	if ($onde(k)>0)
	else
		$onde(k)=0
	endif
k+=1
while(k<n)
endmacro

Window panel2() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(7,57,258,764)
	ListBox list0,pos={1,101},size={248,324},proc=ListPROPreturntoelaborator
	ListBox list0,help={"Display the results of decomposition by increasing delta ppm"}
	ListBox list0,fSize=10,frame=2,listWave=root:list_prop,mode= 1,selRow= 1
	ListBox list0,editStyle= 1,special= {0,15,0},widths={115,49},userColumnResize= 1
	ListBox list1,pos={1,450},size={100,230}
	ListBox list1,help={"Displays the list of patterns whose mass is used to perform linear combination to match the current m/z"}
	ListBox list1,labelBack=(65535,65535,65535),fSize=10,frame=2
	ListBox list1,listWave=root:list_molperm,mode= 1,selRow= 2,special= {0,15,0}
	SetVariable setvar0,pos={1,1},size={250,16},proc=SetVarProcInco,title="m/z to analyse"
	SetVariable setvar0,help={"Displays the current mass over charge used for the calculation belox, step is set to 0.1. It fits to the closest datapoint automatically"}
	SetVariable setvar0,format="%15.20f",limits={-inf,inf,0.1},value= inco,live= 1
	Button button1,pos={163,41},size={80,20},proc=Buttongotopreviousintensity,title="Next taller"
	Button button1,help={"Sets the current m/z to the datapoint immediatly taller than tallest datapoint with same nominal mass as current m/z"}
	Button button1,fColor=(32768,65280,0)
	Button button2,pos={163,61},size={80,20},proc=Buttongotonextintensity,title="Next shorter"
	Button button2,help={"Sets the current m/z to the datapoint immediatly shorter than tallest datapoint with same nominal mass as current m/z"}
	Button button2,fColor=(32768,65280,0)
	Button button3,pos={1,683},size={100,20},proc=ButtonRemoveFromPerm,title="Remove"
	Button button3,help={"Removes selected pattern from list"}
	Button button3,fColor=(65280,43520,0)
	SetVariable setvar1,pos={18,82},size={143,16},bodyWidth=60,title="Displayed results"
	SetVariable setvar1,help={"Sets the length of the proposed formulas list. Please, keep it positive."}
	SetVariable setvar1,limits={0,inf,1},value= nbprop
	ListBox list2,pos={100,450},size={150,230}
	ListBox list2,help={"Sets manual constrains to force or prevent stoechiometic abundance in proposed formulas"}
	ListBox list2,fSize=10,frame=2,listWave=root:minmax,selWave=root:driveminmax
	ListBox list2,mode= 5,special= {0,15,0}
	CheckBox check0,pos={150,683},size={85,14},title="Use Auto Max"
	CheckBox check0,help={"Ignores maximal manual constrains"},value= 1
	Button button4,pos={1,20},size={250,20},proc=ButtonProcheETintense,title="close and tall"
	Button button4,help={"Sets the current m/z above to the tallest datapoint that have the same nominal mass as current m/z."}
	Button button4,fColor=(32768,65280,0)
	Button button5,pos={1,41},size={80,40},proc=Buttoncalculate,title="Calculate"
	Button button5,help={"Calculates combinations of patterns (set below) that could match the current m/z (above) and sort them"}
	Button button5,fColor=(65280,43520,0)
	TitleBox title1,pos={25,435},size={35,13},title="Factors",frame=0
	TitleBox title2,pos={120,435},size={87,13},title="Manual Constrains",frame=0
	Button button0,pos={82,41},size={80,40},proc=ButtonAttribute,title="Attribute"
	Button button0,help={"Opens a panel to run automatic attributions with current settings."}
	Button button0,fColor=(65280,0,0)
EndMacro

Function ButtonAddtoPerm(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			addmoltoperm()
			break
	endswitch

	return 0
End

Function ButtonRemoveFromPerm(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=attributeur list1
			wave/Z list_molperm
			if(numpnts(list_molperm)>0)
			killthismol(v_value)
			calcperm()
			else
			print "Y'en a plus !"
			endif
			break
	endswitch

	return 0
End

function analysemass(bag)
string bag
nvar charge, inco, nbprop, automax
variable masstoanal
wave/Z perm, maxperm
wave/T minmax
razprop()
duplicate/O perm, masstosoustract
masstosoustract=str2num(minmax[x][0])*perm[x]
masstoanal= inco*max(1,abs(charge)) + charge*emass -sum(masstosoustract)
if(itemsinlist(wavelist("molperm_*",";","")) != 0)
	if(automax==1)
		permmax(masstoanal)
		maxperm-=str2num(minmax[x][0])
	else
		permmax(masstoanal)
		maxperm=min(maxperm[x],str2num(minmax[x][1])-str2num(minmax[x][0]))
	endif
endif
trouvemass(masstoanal)
generelistprop(nbprop)
end

////////////
// rework de l'attribution 14 fev 2019
//////
function/wave Bag2perm(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=numpnts(chargebag)
Make/O/D/N=(n)/FREE newperm,index
newperm=target2theo(ciblebag[x][1],ciblebag[x][3])
Makeindex newperm,index
Rotate -1, index
Make/O/D/N=(n,2)/FREE nperm
nperm[][0]=newperm[index[x]]
nperm[][1]=index[x]
return nperm
end

function benchb2p(bag,n)
string bag
variable n
variable t=ticks,k
for(k=0;k<n;k+=1)
	duplicate/O Bag2perm(bag) destroyme
endfor
return (ticks-t)/60.15
end

Function ButtonProcheETintense(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar inco
wave/Z wave1stripped,wave0stripped, incox, incoxdm,incoy
	switch( ba.eventCode )
		case 2: // mouse up
			inco=procheetintense(inco)
			incox=wave0stripped(prochepic(inco,0))
			incoxdm=defmass(incox)
			incoy=wave1stripped(prochepic(inco,0))
			setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
			setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
			break
	endswitch

	return 0
End

function gotoPREVIOUSintensity(pos)//entree en masse reponse en masse
variable pos
wave/Z wave0stripped, dataindex
nvar charge
pos=wave0stripped[prochepic(pos,0)]
findvalue/V=(pos) wave0stripped
findvalue/V=(v_value) dataindex
//pos=round(charge*pos)/charge
//duplicate/O dataindex transclass
//transclass=abs(round(charge*wave0stripped(dataindex))/charge - pos)
//findlevel/Q transclass 0
//killwaves/Z transclass
return wave0stripped(dataindex[v_value-1])
end

function gotoNEXTintensity(pos)
variable pos
wave/Z wave0stripped, dataindex
nvar charge
pos=wave0stripped[prochepic(pos,0)]
findvalue/V=(pos) wave0stripped
findvalue/V=(v_value) dataindex
//pos=round(charge*pos)/charge
//duplicate/O dataindex transclass
//transclass=abs(round(charge*wave0stripped(dataindex))/charge - pos)
//findlevel/Q transclass 0
//killwaves/Z transclass
return wave0stripped(dataindex[v_value+1])
end

Function Buttongotopreviousintensity(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar inco
wave/Z wave1stripped,wave0stripped, incox, incoy
	switch( ba.eventCode )
		case 2: // mouse up
			inco=gotopreviousintensity(inco)
			incox=wave0stripped(prochepic(inco,0))
			incoy=wave1stripped(prochepic(inco,0))
			setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
			setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
			break
	endswitch
	return 0
End

Function Buttongotonextintensity(ba) : ButtonControl
	STRUCT WMButtonAction &ba	
nvar inco
wave/Z wave1stripped,wave0stripped, incox, incoxdm, incoy
	switch( ba.eventCode )
		case 2: // mouse up
			inco=gotonextintensity(inco)
			incox=wave0stripped(prochepic(inco,0))
			incoxdm=defmass(incox)
			incoy=wave1stripped(prochepic(inco,0))
			setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
			setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
			break
	endswitch

	return 0
End





function fullautosimple(nprochains,seuil)// ATTENTION POSSIBLES PROBLEMES DU AU FAIT QUE NBPROP PEUT ETRE SUPERIEUR A J, CONTOURNEMENT PAR STR2NUM !!!
variable nprochains,seuil
variable k, n=nprochains, j 
nvar inco, nbprop
wave/Z dataindex, wave0stripped, deltam, deltai, molecule, transindex, simu_proba, simu_mass, deltappm, test1, labelpro, incox, incoxdm,incoy, wave1stripped, testlocal
wave/T list_prop
string lawave
	for(k=0;k<n;k+=1)
		analysemass("arg")
		Make/O/D/N=(nbprop) critere
		critere=min(seuil, str2num(list_prop[x][1])>seuil)
			if(wavemin(critere)<seuil)
				findvalue/V=(wavemin(critere)) critere
				molecule=0
				tuelesspec()
				lawave="molprop_"+num2str(v_value)
				duplicate/O $lawave Molecule
				genestringmol()
				addthissimu()
			endif
		doupdate
		inco=gotonextintensity(inco)
		incox=wave0stripped(prochepic(inco,0))
		incoxdm=defmass(incox)
		incoy=wave1stripped(prochepic(inco,0))
	endfor
	inco=gotopreviousintensity(inco)
	incox=wave0stripped(prochepic(inco,0))
	incoxdm=defmass(incox)
	incoy=wave1stripped(prochepic(inco,0))
end



function colorit()
variable k, n=itemsinlist(wavelist("molecule_*",";","")),intcum=0, monoch, dich
string lamol, lafor, blaz, ledef
wave/Z trans, index_formules
wave/T list_formules
duplicate/O index_formules color statmag statsou statox statsod
	for(k=1;k<n+1;k+=1)
		lamol="molecule_"+num2str(k)
		lafor="formule_"+num2str(k)
		blaz=list_formules[k-1]
		ledef="DefmassForm_"+num2str(k)
		duplicate/O $lamol trans
		wave laproba=$lafor
		laproba=1e21*exp(-(trans[5]-21)^2/(21))+exp(-(trans[0]-trans[5]*2-29)^2/(29))
		modifygraph/W=agregateur rgb($lafor)=(0,0,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
		modifygraph/W=dmvm rgb($ledef)=(0,0,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
		if (trans[5]!=0)
			modifygraph/W=agregateur rgb($lafor)=(65280,0,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(65280,0,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			color[k-1]=0
			monoch+=1
		endif
		statmag[k-1]=trans[11]
		statox[k-1]=trans[7]
		statsou[k-1]=trans[15]
		statsod[k-1]=trans[10]
		if(sum($lafor)>0)
		intcum+=sum($lafor)
		endif
	endfor
	print monoch, dich
	return intcum
end

function coloritorga()
variable k, n=itemsinlist(wavelist("molecule_*",";","")),intcum=0, monoch, dich
string lamol, lafor, blaz, ledef, lamass
wave/Z trans, index_formules
wave/T list_formules
make/O/N=7 nomb=0, intnomb=0
make/O/T types
types={"CHNOS","CHNO","CHOS","CHNS","CHS","CHO","CHN"}
duplicate/O index_formules color statmag statsou statox statsod
	for(k=1;k<n+1;k+=1)
		lamol="molecule_"+num2str(k)
		lafor="formule_"+num2str(k)
		lamass=lafor+"mass"
		blaz=list_formules[k-1]
		ledef="DefmassForm_"+num2str(k)
		duplicate/O $lamol trans
		if (trans[7]>1 && trans[6]>0)
			print wavemin($lamass)
		endif
		if (trans[0]!=0 && trans[5]!=0 && trans[6]!=0 && trans[7]!=0 && trans[15]!=0)//CHNOS
			modifygraph/W=agregateur rgb($lafor)=(0,0,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(0,0,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,0,0)
			nomb[0]+=1
			intnomb[0]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]!=0 && trans[7]!=0 && trans[15]==0)//CHNO
			modifygraph/W=agregateur rgb($lafor)=(65280,0,52224), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(65280,0,52224), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[1]+=1
			intnomb[1]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]==0 && trans[7]!=0 && trans[15]!=0)//CHOS
			modifygraph/W=agregateur rgb($lafor)=(65280,43520,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(65280,43520,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[2]+=1
			intnomb[2]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]!=0 && trans[7]==0 && trans[15]!=0)//CHNS
			modifygraph/W=agregateur rgb($lafor)=(0,65280,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(0,65280,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[3]+=1
			intnomb[3]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]==0 && trans[7]==0 && trans[15]!=0)//CHS
			modifygraph/W=agregateur rgb($lafor)=(65280,65280,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(65280,65280,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[4]+=1
			intnomb[4]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]==0 && trans[7]!=0 && trans[15]==0)//CHO
			modifygraph/W=agregateur rgb($lafor)=(65280,0,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(65280,0,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[5]+=1
			intnomb[5]+=sum($lafor)
		elseif(trans[0]!=0 && trans[5]!=0 && trans[6]!=0 && trans[7]==0 && trans[15]==0)//CHN
			modifygraph/W=agregateur rgb($lafor)=(0,0,65280), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(0,0,65280), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[6]+=1
			intnomb[6]+=sum($lafor)
		else
			modifygraph/W=agregateur rgb($lafor)=(0,0,0), mode($lafor)=8, Marker($lafor)=23, hideTrace($lafor)=0
			modifygraph/W=dmvm rgb($ledef)=(0,0,0), mode($ledef)=3, Marker($ledef)=18, msize($ledef)=4, useMrkStrokeRGB($ledef)=1,mrkStrokeRGB($ledef)=(65535,65535,65535)
			nomb[7]+=1
			intnomb[7]+=sum($lafor)
		endif
		if (str2num(blaz[strlen(blaz)-4])==1 && sum(trans)!=0)
			monoch+=1
		elseif(str2num(blaz[strlen(blaz)-4])==2 && sum(trans)!=0)
			dich+=1
		endif
		if(sum($lafor)>0)
		intcum+=sum($lafor)
		endif
	endfor
	print monoch, dich
	return intcum
end

function formulebulk()
variable k, n=itemsinlist(wavelist("molecule_*",";","")),intcum=0
string lamol, lafor
wave/T list_formules
Make/O/N=114 moleculebulk=0
	for(k=1;k<n+1;k+=1)
		lamol="molecule_"+num2str(k)
		lafor="formule_"+num2str(k)
		duplicate/O $lamol trans
		moleculebulk+=trans*sum($lafor)		
		intcum+=sum($lafor)
	endfor
	return intcum
end

function bulknonanorh()
variable k, n=itemsinlist(wavelist("molecule_*",";","")),intcum=0
string lamol, lafor
wave/T list_formules
Make/O/N=114 moleculebulk=0
	for(k=1;k<n+1;k+=1)
		lamol="molecule_"+num2str(k)
		lafor="formule_"+num2str(k)
		duplicate/O $lamol trans
		if (trans[10]==1)
			trans[10]+=-1
		else
			trans[0]+=-1
		endif
		moleculebulk+=trans*sum($lafor)		
		intcum+=sum($lafor)
	endfor
	return intcum
end

function spotbadpic()
wave/Z bad_peak, index_formules
variable k,n = numpnts(bad_peak)
string lafor
for(k=0;k<n;k+=1)
findvalue/V=(bad_peak[k])/T=0.001 index_formules
lafor="formule_"+num2str(v_value+1)
print bad_peak[k],lafor
modifygraph/W=agregateur mode($lafor)=8, opaque($lafor)=1
endfor
end

function spotnitrotoox()
wave/Z CO5H6asN7, index_formules
variable k,n = numpnts(CO5H6asN7)
string lafor
for(k=0;k<n;k+=1)
findvalue/V=(CO5H6asN7[k])/T=0.001 index_formules
lafor="formule_"+num2str(v_value+1)
print CO5H6asN7[k],lafor
modifygraph/W=agregateur mode($lafor)=8, opaque($lafor)=1, marker($lafor)=23
endfor
end

function spotabsent()
wave/Z absent, index_formules
variable k,n = numpnts(absent)
string lafor
for(k=0;k<n;k+=1)
findvalue/V=(absent[k])/T=0.001 index_formules
lafor="formule_"+num2str(v_value+1)
print absent[k],lafor
modifygraph/W=agregateur mode($lafor)=8, opaque($lafor)=1, marker($lafor)=8
endfor
end

function topbright()
variable binf, bsup
wave/Z topbrighty, topbrightx
getaxis/W=agregateur/Q bottom
binf=prochepic(v_min,0)-1
bsup=prochepic(v_max,0)+1
duplicate/O wave1stripped transindex
duplicate/O wave0stripped transx
deletepoints 0, binf, transindex, transx
deletepoints abs(bsup-binf), numpnts(transindex), transindex, transx
sort/R transindex, transindex, transx
topbrighty=transindex[x]
topbrightx=transx[x]
killwaves transindex, transx
end

function recal()
string lesmasses=wavelist("*_0stripped",";",""), lamasse
wave/Z trans
variable k, n=itemsinlist(lesmasses), poscaf, shift
	for(k=0;k<n;k+=1)
		lamasse=stringfromlist(k,lesmasses)
		duplicate/O $lamasse trans
		findlevel/Q trans,195.0876552
		poscaf=trans[round(v_levelx)]
		shift=poscaf-195.0876552
		trans-=shift
		duplicate/O trans $lamasse
		printf "%+15.7f\r", trans[round(v_levelx)]
	endfor
end

macro defmassplot(defmass,mass,color)
string defmass,mass,color
	PauseUpdate; Silent 1		// building window...
	Display /W=(85.5,306.5,538.5,666.5) $defmass vs $mass
	ModifyGraph margin(right)=100,width=300,height=300,wbRGB=(65280,0,26112),gbRGB=(0,0,0)
	ModifyGraph mode=2
	ModifyGraph zColor($defmass)={$color,*,*,BlackBody}
	ModifyGraph mirror=1
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(left)=3
	ModifyGraph axOffset(left)=-1
	ModifyGraph axThick=2
	ModifyGraph axRGB=(32768,54528,65280)
	ModifyGraph tlblRGB=(32768,54528,65280)
	ModifyGraph alblRGB=(32768,54528,65280)
	ModifyGraph lblLatPos(left)=-2
	Label left "\\f03Mass defect (u)"
	Label bottom "m/z"
	ColorScale/C/N=text0/G=(32768,54528,65280)/B=(65280,0,26112)/A=MC/X=67.75/Y=0.50
	ColorScale/C/N=text0 trace=$defmass, fstyle=1, lblRot=180
	AppendText "\\f03Relative intensity (log scale)"
	ModifyGraph fSize=20
	ColorScale/C/N=text0 "\\Z20\\f03Relative intensity (log scale)"
	SetAxis left -0.5,0.5
	setaxis bottom 50, 500
EndMacro
Window Table0() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,41.75,510,235.25) defmassreliso12,defmassreliso23,defmassreliso13
	ModifyTable format(Point)=1
EndMacro

Window dmvm() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(951,42.5,1389,539.75) simu_probadm vs simu_mass
	AppendToGraph roiy vs roix
	AppendToGraph yplode vs xplode
	AppendToGraph cauchyyplode vs cauchyxplode
	AppendToGraph currentseg[*][0] vs currentseg[*][1]
	AppendToGraph wave2stripped vs wave0stripped
	AppendToGraph roi2 vs roi0
	AppendToGraph defmasstree vs tree
	AppendToGraph incoxdm vs incox
	ModifyGraph gbRGB=(47872,47872,47872)
	ModifyGraph mode(simu_probadm)=3,mode(wave2stripped)=3,mode(roi2)=3,mode(defmasstree)=3
	ModifyGraph mode(incoxdm)=3
	ModifyGraph marker(simu_probadm)=18,marker(wave2stripped)=19,marker(roi2)=19,marker(defmasstree)=55
	ModifyGraph marker(incoxdm)=42
	ModifyGraph lSize(yplode)=1.2,lSize(cauchyyplode)=1.2,lSize(currentseg)=2,lSize(roi2)=1.2
	ModifyGraph rgb(simu_probadm)=(0,52224,52224),rgb(yplode)=(0,0,65280),rgb(currentseg)=(0,52224,52224)
	ModifyGraph rgb(wave2stripped)=(0,0,0),rgb(roi2)=(36864,14592,58880),rgb(defmasstree)=(34816,34816,34816)
	ModifyGraph rgb(incoxdm)=(26368,0,52224)
	ModifyGraph msize(wave2stripped)=5,msize(roi2)=4,msize(defmasstree)=7,msize(incoxdm)=5
	ModifyGraph mrkThick(wave2stripped)=0.001,mrkThick(defmasstree)=1
	ModifyGraph opaque(wave2stripped)=1,opaque(incoxdm)=1
	ModifyGraph hideTrace(yplode)=1,hideTrace(cauchyyplode)=1
	ModifyGraph useMrkStrokeRGB(simu_probadm)=1,useMrkStrokeRGB(wave2stripped)=1,useMrkStrokeRGB(roi2)=1
	ModifyGraph useMrkStrokeRGB(defmasstree)=1
	ModifyGraph mrkStrokeRGB(simu_probadm)=(65535,65535,65535)
	ModifyGraph zmrkSize(wave2stripped)={logwave1,*,*,1,5},zmrkSize(defmasstree)={treeerror,*,10,6,0}
	ModifyGraph zColor(wave2stripped)={logwave1,-0.0772956579325359,*,YellowHot},zColor(defmasstree)={treeColor,*,*,Web216,1}
	ModifyGraph zColorMin(wave2stripped)=(65535,65535,65535)
	ModifyGraph zero(left)=1
	ModifyGraph mirror=1
	ModifyGraph fStyle=1
	ModifyGraph axOffset(left)=-2.14286,axOffset(bottom)=-0.388889
	Label left "Mass defect"
	Label bottom "m/z"
	SetAxis left -0.499888000000055,0.499706999999944
	SetAxis bottom 150.002339,699.521952
	ControlBar 40
	ValDisplay valdisp0,pos={6,1},size={89,14},title="ROI size"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"numpnts(roi1)"
	ValDisplay valdisp1,pos={106,1},size={99,14},title="IntCum %"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"round(100*(sum(roi1)/sum(wave1stripped)))"
	ValDisplay valdisp2,pos={326,1},size={74,14},title="InPnts"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000},value= #"twosig"
	SetVariable setvar0,pos={205,1},size={120,16},proc=Proc_setintenscumthres,title="IntCumThres"
	SetVariable setvar0,limits={0,1,0.1},value= intenscumthres
	SetVariable setvar1,pos={6,18},size={113,16},proc=SetVarProc_6,title="Depth of tree"
	SetVariable setvar1,limits={0,inf,1},value= _NUM:0
	ValDisplay valdisp3,pos={124,18},size={142,14},title="Segment length :"
	ValDisplay valdisp3,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp3,value= #"sqrt((currentseg[0][1]-currentseg[1][1])^2+(currentseg[0][0]-currentseg[1][0])^2)"
	PopupMenu popup0,pos={269,18},size={123,21},proc=PopMenuProc_3,title="Use this list:"
	PopupMenu popup0,mode=25,popvalue="Factors",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	SetWindow kwTopWin,hook=getdmvmMSMSmode,hookevents=3
	SetWindow kwTopWin,hook(dmvm)=getdmvmMSMSmode
EndMacro



function plotles(lesy,lesx)
string lesy, lesx
display
variable k, n=itemsinlist(wavelist(lesy,";",""))
string ondey, ondex
	if (stringmatch(lesx,"")==1)
		for(k=0;k<n;k+=1)
		ondey=stringfromlist(k,wavelist(lesy,";",""))
		appendtograph $ondey
	endfor
	else
		for(k=0;k<n;k+=1)
	ondey=stringfromlist(k,wavelist(lesy,";",""))
	ondex=stringfromlist(k,wavelist(lesx,";",""))
		appendtograph $ondey vs $ondex
		print numpnts($ondey)
	endfor
	setaxis bottom wavemax($ondex), wavemin($ondex)
	endif
end

Function SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
nvar top2dis
	switch( sva.eventCode )
		case 1: // mouse up
		 topbright()
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			make/O/N=(top2dis) topbrighty, topbrightx
			topbright()
			break
	endswitch

	return 0
End


Function agregaReset(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			setaxis/W=agregateur/A bottom// wavemin(wave0stripped), wavemax(wave0stripped)
			setaxis/W=dmvm/A bottom
			DoUpdate
			setaxis/W=agregateur /A=2 left
			topbright()
			calcinfo()
			break
	endswitch

	return 0
End

Function ButtonAttribute(ba) : ButtonControl
	STRUCT WMButtonAction &ba
string lespanel3=winlist("panel3*",";",""), lepanel3
variable k, n=itemsinlist(lespanel3)
	switch( ba.eventCode )
		case 2: // mouse up
			if(stringmatch(lespanel3,""))
			Execute "panel3()"
			else
				for(k=0;k<n;k+=1)
					lepanel3=stringfromlist(k,lespanel3)
					killwindow $lepanel3
				endfor
			endif
			break
	endswitch

	return 0
End

Function ButtonGoForAttribution(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar nexttoatt, ppmthreshold
	switch( ba.eventCode )
		case 2: // mouse up
				fullauto(nexttoatt,ppmthreshold)
			break
	endswitch

	return 0
End

Window showprogbar() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(237,349.25,733.5,398)
	AppendImage/T progressbarre
	SetAxis/R left 1,-0.5
	ModifyImage progressbarre ctab= {*,*,Grays,0}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph tick=3
	ModifyGraph mirror=2
	ModifyGraph nticks(top)=3
	ModifyGraph minor=1
	ModifyGraph noLabel=2
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	SetAxis/A/R left
EndMacro

function findclosepntagreg(xx,yy)
variable xx,yy
nvar inco
wave/Z Wave0stripped, wave1stripped, incox,incoy,incoxdm
variable stretchx, stretchy
duplicate/O wave0stripped distance
getaxis/Q/W=agregateur bottom
stretchx=v_max-v_min
getaxis/Q/W=agregateur left
stretchy=v_max-v_min
distance=( (wave0stripped-xx)/stretchx )^2 + ( (wave1stripped-yy)/stretchy )^2
findvalue/V=(wavemin(distance)) distance
incox=wave0stripped[v_value]
incoy=wave1stripped[v_value]
incoxdm=defmass(incox)
inco=incox[0]
killwaves distance
end

function findclosepntdmvm(xx,yy)
variable xx,yy
nvar inco
wave/Z Wave0stripped, wave2stripped, incox,incoy,incoxdm, wave1stripped
variable stretchx, stretchy
duplicate/O wave0stripped distance
getaxis/Q/W=dmvm bottom
stretchx=v_max-v_min
getaxis/Q/W=dmvm left
stretchy=v_max-v_min
distance=( (wave0stripped-xx)/stretchx )^2 + ( (wave2stripped-yy)/stretchy )^2
findlevel/Q distance (wavemin(distance))
incox=wave0stripped[v_levelx]
incoy=wave1stripped[v_levelx]
incoxdm=defmass(incox)
inco=incox[0]
killwaves distance
end

function loadbasic()
string name="Sample#_"+num2str(abs(enoise(1e12)))
prompt name, "Name this Sample"
doprompt "Enter a good name for this sample",name
if(v_flag==0)
	string cheminversfichier
	string nomfichier, fichiery, fichierx
	GetFileFolderInfo/Q
	NewPath/O/Q ledossier ParseFilePath(1, S_Path, ":", 1, 0)
	nomfichier=ParseFilePath(0, S_Path, ":", 1, 0)
	LoadWave/O/A/Q/G/D/P=ledossier nomfichier
	nomfichier=name[0,19]
	fichiery=nomfichier+"_1stripped"
	fichierx=nomfichier+"_0stripped"
	Rename wave0, $fichierx
	Rename wave1, $fichiery
endif
end

Function ButtonLoader(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			loadbasic()
			listelesdonnees()
			break
	endswitch

	return 0
End

Function ButtonChromaLoader(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			loadchroma()
			break
	endswitch
	return 0
End

Function ButtonChromaKill(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			killlechroma()
			break
	endswitch
	return 0
End

Function ButtonKill(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			killladata()
			break
	endswitch
	return 0
End

function killladata()
nvar inco
controlinfo/W=advancedmanager TabCon_List_DATA
wave/T listwave=$S_value
string laun,lazero, lehistoiso, lehistoslo, linfo
variable taille
wave/Z liste_donnees
if(numpnts(liste_donnees)>0)
laun=listwave[V_value][1]+"_1stripped"
lazero=listwave[V_value][1]+"_0stripped"
lehistoiso=listwave[V_value][1]
linfo=lehistoiso+"_info"
lehistoiso="Histoiso_"+lehistoiso
//print lehistoiso
lehistoslo=listwave[V_value][1]
lehistoslo="Histoslo_"+lehistoslo
//print lehistoslo
killwaves/Z $laun, $lazero,$lehistoslo,$lehistoiso, $linfo
listelesdonnees()
endif
end

function killlechroma()
nvar inco
controlinfo/W=advancedmanager TabCon_List_CHROMA
wave/T listwave=$S_value
string chrom2kill=wavelist(listwave[V_value][1]+"_*chroma",";",""),target
if(numpnts(liste_chroma)>0)
variable k, n=itemsinlist(chrom2kill)
for(k=0;k<n;k+=1)
	target=stringfromlist(k,chrom2kill)
	killwaves/Z $target
endfor
listelesdonnees()
endif
end


Function ButtonDupeRname(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=advancedmanager TabCon_List_DATA
			wave/T listwave=$S_value
			string newname= listwave[V_value][1]+"_"
			prompt newname, "New Name"
			Doprompt "Enter a new good name for this sample",newname
			newname=UniqueName(newname[0,18], 1, 0)
			if(v_flag==0)
				renomeactive(newname)
			endif
			break
	endswitch

	return 0
End


function recalonthefly(intstd)
variable intstd
variable posintstd, shift
nvar inco
wave/Z wave0stripped
findlevel/Q wave0stripped, intstd
posintstd=wave0stripped[round(v_levelx)]
shift=posintstd-intstd
wave0stripped-=shift
end 

function hardrecal(intstd)
variable intstd
variable shift
nvar inco
wave/Z wave0stripped
shift=inco-intstd
wave0stripped-=shift
end 

Function ButtonCalib(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z simu_mass, simu_proba
			variable xeumap=wavemax(simu_proba), mass2recal
			findvalue/V=(xeumap) simu_proba
			mass2recal=simu_mass[v_value]
			recalonthefly(mass2recal)
			break
	endswitch

	return 0
End

Function ButtonHardCalib(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z simu_mass, simu_proba
			variable xeumap=wavemax(simu_proba), mass2recal
			findvalue/V=(xeumap) simu_proba
			mass2recal=simu_mass[v_value]
			hardrecal(mass2recal)
			break
	endswitch

	return 0
End

function fillvol(kprem)
variable kprem
wave/Z M_3DVertexList, w_cross, wave0stripped, wave1stripped, wave2stripped, dataindex
make/O/N=(kprem,3) fortrig
fortrig[][0]=wave0stripped[dataindex[x]]
fortrig[][1]=0
fortrig[][2]=wave2stripped[dataindex[x]]
insertpoints 0,1, fortrig
fortrig[0][1]=1
triangulate3D fortrig
//killwaves fortrig, vecta,vectb,vectc,vectd, volumes, cote
end

function sortlines(mat)
wave/Z mat
variable k, n=dimsize(mat,0), m=dimsize(mat,1)
make/O/N=(m) trans
	for(k=0;k<n;k+=1)
		trans=mat[k][x]
		sort/R trans trans
		mat[k][]=trans[y]
	endfor
end

function plotdelau(mat,mat2)
wave/Z mat, mat2
make/O/N=0 xdelau, ydelau, segments, pentes
variable k, n=dimsize(mat,0), m=dimsize(mat,1)-1,i
	for(k=0;k<n;k+=1)
		for(i=0;i<m;i+=1)
			insertpoints 0, 1, xdelau, ydelau
			xdelau[0]=mat2[mat[k][i]][0]
			ydelau[0]=mat2[mat[k][i]][2]
		endfor
		insertpoints 0,3,segments, pentes
		segments[0]=sqrt((xdelau[1]-xdelau[0])^2+(ydelau[1]-ydelau[0])^2)
		segments[1]=sqrt((xdelau[2]-xdelau[1])^2+(ydelau[2]-ydelau[1])^2)
		segments[2]=sqrt((xdelau[2]-xdelau[0])^2+(ydelau[2]-ydelau[0])^2)
		pentes[0]=(xdelau[1]-xdelau[0])/(ydelau[1]-ydelau[0])
		pentes[1]=(xdelau[2]-xdelau[1])/(ydelau[2]-ydelau[1])
		pentes[2]=(xdelau[2]-xdelau[0])/(ydelau[2]-ydelau[0])
		insertpoints 0,1,xdelau, ydelau
		xdelau[0]=NaN
	endfor
end

function choperoi()
wave/Z wave0stripped, wave2stripped, roix, roiy,w_findlevels, wave1stripped
variable k, n=numpnts(wave0stripped), j
make/O/N=0 roipnts
if(numpnts(roiy)!=0)
	for(k=0;k<n;k+=1)
		findlevels/Q roiy wave2stripped[k]
		duplicate/O w_findlevels paroix
		paroix=roix[w_findlevels[x]]-wave0stripped[k]<0
		findlevels/Q roix wave0stripped[k]
		duplicate/O w_findlevels paroiy
		paroiy=roiy[w_findlevels[x]]-wave2stripped[k]<0
		//print mod(sum(paroix),2),mod(sum(paroiy),2)
		if(mod(sum(paroix),2)!=0 && mod(sum(paroiy),2)!=0)
			insertpoints 0,1, roipnts
			roipnts[0]=k
		endif
	endfor
	duplicate/O roipnts roi1 roi0 roi2
	roi1=wave1stripped[roipnts[x]]
	roi0=wave0stripped[roipnts[x]]
	roi2=wave2stripped[roipnts[x]]
else
make/O/N=0 roi1, roi0,roi2
make/O/N=0 roipnts
endif
calcinfo()
end

function simupente(xeuma)
variable xeuma
wave/Z perm, wave0stripped, wave1stripped
variable n=numpnts(perm), k, m, maxx, maxy,contrainte
duplicate/O perm maxpermpente prodpente
maxpermpente=2*ceil(xeuma/perm)+1
if(xeuma<wavemin(perm))//autotune l'inflation en cas de xeuma plus petit que les factors
	duplicate/O perm temperm
	sort perm temperm
	differentiate/METH=1 temperm
	temperm=ceil(abs((xeuma/temperm)/2))
	Rotate -1, temperm
	maxpermpente+=(2*(temperm))
	killwaves temperm
endif
prodpente=1
for(k=0;k<n-1;k+=1)
	prodpente[k+1]=prodpente[k]*maxpermpente[k]
endfor
m=prodpente[numpnts(prodpente)-1]*maxpermpente[numpnts(prodpente)-1]
make/O/N=(m,n) criblepente
criblepente = (mod(trunc(x/prodpente[y]),maxpermpente[y]))-((maxpermpente[y]-1)/2)
make/O/N=(m) pentes_simu, modules_simu
pentes_simu=0
for(k=0;k<n;k+=1)
	pentes_simu+=criblepente[x][k]*perm[k]//calcul la masse virtuelle (à indices neg)
endfor
modules_simu=sqrt((pentes_simu)^2+(pentes_simu-round(pentes_simu))^2)
pentes_simu=1-round(pentes_simu)/pentes_simu
//make/O histoslo_simu
//histogram/B={-0.01,0.02/1000,1000} pentes_simu histoslo_simu
//killwaves pentes_simu, criblepente, maxpermpente, prodpente
//duplicate/O histoslo_simu xplode_simu yplode_simu
//xplode_simu=(histoslo_simu)*cos(atan(pnt2x(histoslo_simu,p)))
//yplode_simu=(histoslo_simu)*sin(atan(pnt2x(histoslo_simu,p)))
//maxx=max(wavemax(xplode_simu), abs(wavemin(xplode_simu)))
//maxy=max(wavemax(yplode_simu), abs(wavemin(yplode_simu)))
//contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
//xplode_simu*=contrainte
//yplode_simu*=contrainte
end

function pentestocha(n,m)
variable n,m
wave/Z perm
variable k=0,maxx,maxy, contrainte
duplicate/O perm combistock maxpermstocha
maxpermstocha=ceil(m/perm)+1
make/O/N=1000 histo_trans, histoslo_stocha
make/O/N=1 pente_stocha
histogram/B={-0.01,0.02/1000,1000} pente_stocha histoslo_stocha
histoslo_stocha=0
do
combistock=round(enoise(maxpermstocha))*perm
if(sum(combistock)<m)
pente_stocha=1-round(sum(combistock))/(sum(combistock))
histogram/B={-0.01,0.02/1000,1000} pente_stocha histo_trans
histoslo_stocha+=histo_trans
//doupdate
k+=1
endif
while(k<n)
duplicate/O histoslo_stocha xplode_stocha yplode_stocha
xplode_stocha=(histoslo_stocha)*cos(atan(pnt2x(histoslo_stocha,p)))
yplode_stocha=(histoslo_stocha)*sin(atan(pnt2x(histoslo_stocha,p)))
maxx=max(wavemax(xplode_stocha), abs(wavemin(xplode_stocha)))
maxy=max(wavemax(yplode_stocha), abs(wavemin(yplode_stocha)))
contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
xplode_stocha*=contrainte
yplode_stocha*=contrainte
end


Function ButtonProc_1(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			SetAxis bottom 1.995,2.005
			SetAxis top 1.995,2.005
			SetAxis top2 1.995,2.005
			SetAxis top3 1.995,2.005
			break
	endswitch

	return 0
End

Function ButtonZ2(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			SetAxis bottom 1.995,2.005
			SetAxis top 1.995,2.005
			SetAxis top2 1.995,2.005
			SetAxis top3 1.995,2.005
			break
	endswitch

	return 0
End

Function ButtonZ1(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			SetAxis bottom 0.995,1.005
			SetAxis top 0.995,1.005
			SetAxis top2 0.995,1.005
			SetAxis top3 0.995,1.005
			break
	endswitch

	return 0
End

Function ButtonZ0(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			setaxis/A
			break
	endswitch

	return 0
End

function fakeHCO(m)
variable m
wave/Z Example_Brut_0stripped, mendev_masses
duplicate/O Example_Brut_0stripped h
h=abs(round(enoise(round(m/14))))
Example_Brut_0stripped=h*(mendev_masses[5][0]+2*mendev_masses[0][0])+2*(1-min(h,abs(round(gnoise(h)))))*mendev_masses[0][0]+(min(h,abs(round(gnoise(h)))))*mendev_masses[7][0]
end

function fakeHC(m)
variable m
wave/Z Example_Brut_0stripped, mendev_masses
make/O/N=88058 Example_Brut_0stripped, Example_Brut_1stripped, h
h=abs(round(enoise(round(m/14))))
Example_Brut_1stripped=1
Example_Brut_0stripped=h*(mendev_masses[5][0]+2*mendev_masses[0][0])+2*(1-min(h,abs(round(gnoise(h)))))*mendev_masses[0][0]
killwaves h
end

Window Report() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(226,161,1049,755)
	ListBox list0,pos={1,2},size={821,592},listWave=root:Rapport,userColumnResize= 1
EndMacro

function generateSlopeChart()
wave/Z perm, wave0stripped
wave/T list_molperm
variable maxx, k,n=numpnts(perm)
string etik
duplicate/O perm ychart xchart
xchart=sqrt(abs(1/((1*(wavemax(wave0stripped)-wavemin(wave0stripped)))^(-2) + ((1-round(perm)/perm)/0.5)^2) ))
ychart=xchart*(1-round(perm)/perm)
duplicate/O xchart xchartmark
duplicate/O ychart ychartmark
display ychartmark vs xchartmark
setaxis left, -0.5, 0.5
setaxis/A bottom 0,wavemax(wave0stripped)
ModifyGraph mode(ychartmark)=3,marker(ychartmark)=46
	for(k=0;k<n;k+=1)
		insertpoints 3*k, 2, ychart, xchart
		ychart[3*k+1]=0
		xchart[3*k+1]=0
		ychart[3*k]=Nan
		xchart[3*k]=Nan
		etik="text"+num2str(k)
		TextBox/C/N=$etik/M/H=15/A=LC/X=(100*xchartmark[k]/((wavemax(wave0stripped))))/Y=(100*ychartmark[k]) list_molperm[k]
	endfor
appendtograph ychart vs xchart
AppendToGraph rawyplode vs rawxplode
ModifyGraph muloffset(rawyplode)={0.2,0.2}
ModifyGraph rgb(rawyplode)=(0,0,65280)
end

Function TabSettings()
nvar nexttoatt, nbprop, ppmthreshold, charge, krikriter
Wave/T list_molperm, minmax
variable m=numpnts(list_molperm)
Make/O/T/N=(5+m,3) Settings
Settings[0][0]="Next m/z positions to analyse :"
Settings[0][1]=num2str(nexttoatt)
Settings[1][0]="For each, check these first :"
Settings[1][1]=num2str(nbprop)
Settings[2][0]="With this threshold (ppm) :"
Settings[2][1]=num2str(ppmthreshold)
Settings[3][0]="With this ion charge :"
Settings[3][1]=num2str(charge)
Settings[4][0]="And this combo probability :"
Settings[4][1]=num2str(krikriter)
Settings[5,m+5][0]=list_molperm[x-5]
Settings[5,m+5][1,2]=minmax[x-5][y-1]
//Edit Settings
End

function killmissmatch()
variable k, n=itemsinlist(wavelist("molecule_*",";",""))
string lamol, lafor, blaz, ledef, lamas
wave/Z trans, index_formules
wave/T list_formules
duplicate/O index_formules color statmag statsou statox statsod
	for(k=n;k>0;k-=1)
		lamol="molecule_"+num2str(k)
		lamas="formule_"+num2str(k)+"mass"
		lafor="formule_"+num2str(k)
		ledef="DefmassForm_"+num2str(k)
		duplicate/O $lamol trans
		trans=abs(trans)
		if(sum(trans)==0)
		Removefromgraph/W=agregateur $lafor
		Removefromgraph/W=dmvm $ledef
		Killwaves/Z $lafor,$lamas, $lamol, $ledef
		deletepoints/M=0 k-1, 1, list_formules, index_formules, list_targets
		MolbagWithdrawFrom({k-1},"Current")
		endif
	endfor
	majformulewavesnames()
	majliste()
	majlegende()
end

function meltattrib(name,term)
string name,term
string lesx=wavelist("Formule_*mass",";",""), lesy=replacestring("mass",lesx,""), newy, newx
if(strlen(lesx))
	newy=name+"_1"+term
	newx=name+"_0"+term
	concatenate/O/NP lesy, $newy
	concatenate/O/NP lesx, $newx
	sort $newx, $newy, $newx
endif
//listelesdonnees()
end

//===========ARIZONA DREAM STARTS=============

///// getslopestat et generelesmats corrigées pour introduire les centremasses pour recalibration le 7 novembre 2013

function attribCombi(minimod,slo)
variable minimod,slo
variable i,j,fac
Variable/G lamasse
svar nommol
string transname
nvar charge
minimod*=abs(charge)
simupente(minimod)
wave/Z criblepente, pentes_simu, modules_simu, maxpermpente, molecule
modules_simu=abs(modules_simu-minimod)
pentes_simu=abs(pentes_simu-slo)
duplicate/O modules_simu indexmod indexslo
makeindex modules_simu indexmod
makeindex pentes_simu indexslo
duplicate/O molecule molpente
i=numpnts(maxpermpente)
molpente=0
nommol=""
//if(pentes_simu[indexmod[0]]>0 && pentes_simu[indexmod[0]]<2e-5 || slo==pi)
	for(j=0;j<i;j+=1)
		transname="molperm_"+num2str(j)
		Duplicate/O $transname, transmolprop
		molpente+=criblepente[indexmod[0]][j]*transmolprop
	endfor
	duplicate/O molpente molecule
	fac=simplimol()
	molecule/=fac
	lamasse=genestringpente()
if(pentes_simu[indexmod[0]]>=2e-5)
nommol+="   Slo_Miss"
endif
print modules_simu[indexmod[0]], pentes_simu[indexmod[0]]
sort modules_simu, modules_simu, pentes_simu
end

function genestringpente()
wave/Z molecule, nb_isotopes,mendev_masses, mendev_abundances
wave/T elements
variable n=numpnts(elements), k,i, j, nbmass=0, massmol=0, intmass=0, lamasse
string elem
string/G nommol
nommol=""
	for(k=0;k<n;k+=1)//verifie la presence de chaque element
		if (molecule[k]!=0) // si l' element k est present
			elem=elements[k]
			nommol=nommol+elem+"\\B"+num2str(molecule[k])+"\M"//genere la formule brute
		endif
	endfor
duplicate/O molecule molemass
molemass=molecule[x]*mendev_masses[x][0]
lamasse=sum(molemass)
killwaves/Z molemass
return lamasse
end

function findmostzero()
variable para,k,n,a=0,b,t=ticks
wave/Z modtrans, pentrans,interceptrans
duplicate/O modtrans modtransmod rapportmod rapportslo
rapportmod=0
do
differentiate rapportslo
//differentiate pentrans
//rapportslo*=2
duplicate/O rapportslo modtransmod
modtransmod=abs(modtransmod)
//sort modtransmod modtransmod
//simplimod("modtransmod")
//rapportslo=0
//n=numpnts(modtransmod)
//for(k=0;k<n;k+=1)
	para=wavemin(modtransmod)
	rapportmod+=(abs(rapportslo)<=para)
	//doupdate
//endfor
b=sum(rapportmod)-a
a=sum(rapportmod)
doupdate
while(b!=2 && sum(rapportslo)!=0)
duplicate/O rapportmod modtransmod
differentiate modtransmod
findlevels/Q modtransmod 0
duplicate/O w_findlevels priosep priopoint priolevel prioindex priomod prioslo priointercept
priolevel=rapportmod[priopoint]/wavemax(rapportmod)
makeindex/R priolevel prioindex
sort/R priolevel priolevel
priomod=modtrans[priopoint[prioindex]]
priosep=modtrans[priopoint[prioindex+1]]-modtrans[priopoint[prioindex-1]]//priomod
prioslo=pentrans[priopoint[prioindex]]
priointercept=interceptrans[priopoint[prioindex]]
//rapportslo=rapportmod(w_findlevels)*(modtransmod(floor(w_findlevels))>0 && modtransmod(ceil(w_findlevels))<0)
//print ticks-t
end

function simplifysmatch()
wave/Z slopematch, checkslo
variable k,i, n=dimsize(slopematch,0)
duplicate/O checkslo checksloindex
makeindex/R checkslo checksloindex
for(k=0;k<n;k+=1)
	for(i=0;i<n;i+=1)
		if(slopematch[i][checksloindex[k]]==1)
			slopematch[][i]=0
		endif
	endfor
endfor
end

function simplimol()
wave/Z molecule,mendev_masses
variable k,i,n=numpnts(molecule),a,b
duplicate/O molecule molemass
molemass=molecule[x]*mendev_masses[x][0]
make/O/N=(n,n) divcrible
divcrible=1e12
for(k=0;k<n;k+=1)
	for(i=0;i<n;i+=1)
		a=molecule[k]
		b=molecule[i]
		if(a!=0 && b!=0)
			divcrible[k][i]=gcd(a,b)
		endif
	endfor
endfor
a=wavemin(divcrible)*sign(mean(molemass))
killwaves divcrible, molemass
if(a<=wavemax(molecule))
	return sign(a)
else
	return 1
endif
end

function spanslopes(kprem)
variable kprem
variable k,n=1000,t
svar nommol
wave/Z slopetrix
generelesmat(kprem)
make/O/N=(n) histo
histo=0
display histo
newimage slopetrix
wave/Z modtrans,priomod,prioslo, molecule, histo
for(k=0;k<n/2;k+=1)
	t=ticks
	getslopestat(k*2e-5)
	findmostzero()
	analysemult(-1)
	histo[k+n/2]=ticks-t

	t=ticks	
	getslopestat((-k-1)*2e-5)
	findmostzero()
	analysemult(-1)
	histo[n/2-1-k]=ticks-t
	doupdate
endfor
end

function simplimod(lawave)
string lawave
duplicate/O $lawave simplified
variable k,i=0,n=numpnts(simplified)
make/O/N=1 rapportsimp
rapportsimp[0]=simplified[0]
for(k=1;k<n;k+=1)
	if(simplified[k]!=simplified[k-1])
	insertpoints 0,1,rapportsimp
	rapportsimp[0]=simplified[k]
	endif
endfor
sort rapportsimp rapportsimp
duplicate/O rapportsimp $lawave
killwaves rapportsimp simplified
end

function analysemult(kprio)//-1=all
variable kprio
svar nommol
wave/Z priomod, prioslo, priosep, priolevel, molecule
variable k,i,n
if(kprio==-1)
n=numpnts(priomod)
else
n=kprio
endif
make/O/N=0 unikmod, uniksep, uniklvl
simupente(0)
wave/Z modules_simu
for(k=0;k<n && priolevel[k]!=0;k+=1)
	//doupdate
	duplicate/O unikmod uniktest
	uniktest=abs(unikmod-priomod[k])<uniksep+priosep[k]
		if(sum(uniktest)==0 && priomod[k]<2000 || k==0)
			attribcombi(priomod[k],prioslo[k])
				if(stringmatch(nommol,"Slope_missmatch")==0)
				insertpoints 0,1,unikmod,uniksep,uniklvl
				unikmod[0]=priomod[k]
				uniksep[0]=modules_simu[0]//+priosep[k]
				uniklvl[0]=priolevel[k]
				filltrix(unikmod[0],uniksep[0])
				endif
		elseif(sum(uniktest)==0 && priomod[k]>=2000 && k!=0)//PROBLEME ISOBARISME GRRRRR tweak avec 2000
			attribcombi(abs(priomod[k]-unikmod[0]),pi)
				if(stringmatch(nommol,"Slope_missmatch")==0)
				insertpoints 0,1,unikmod,uniksep,uniklvl
				unikmod[0]=priomod[k]
				uniksep[0]=modules_simu[0]//+priosep[k]+uniksep[0]
				uniklvl[0]=priolevel[k]
				filltrix(unikmod[0],uniksep[0])
				endif
		endif
endfor
end

function filltrix(lemod,lasep)
variable lemod,lasep
nvar lamasse
variable k,n,m,lemodcalc
svar nommol
wave/Z slopetrix, modtrans, rapportslo,modtransmod,xxtrans,w_findlevels, molecule,calcerror, modules
wave/T slopenames
rapportslo=abs(modtrans-lemod)<=lasep/2
findlevels/Q rapportslo 1
duplicate/O w_findlevels modtransmod
modtransmod=xxtrans[w_findlevels]
n=numpnts(modtransmod)
m=dimsize(slopetrix,0)
duplicate/O molecule checksum
lemodcalc=sqrt((defmass(lamasse))^2+(lamasse)^2)
for(k=0;k<n;k+=1)
	checksum=slopetrix[mod(modtransmod[k],m)][floor(modtransmod[k]/m)][x]
	if(sum(checksum)!=0)
	//print slopenames[mod(modtransmod[k],m)][floor(modtransmod[k]/m)]
	else
	slopetrix[mod(modtransmod[k],m)][floor(modtransmod[k]/m)][]=molecule[z]
	slopenames[mod(modtransmod[k],m)][floor(modtransmod[k]/m)]=nommol
	calcerror[mod(modtransmod[k],m)][floor(modtransmod[k]/m)]=abs(modules[mod(modtransmod[k],m)][floor(modtransmod[k]/m)]-lemodcalc)
	endif
endfor
end

function libere()
killwaves/Z coefz,cristcheck,slopenames,moletest,checksloindex,histo,calcerror,checksum
killwaves/Z slopetrix, indexperm, test1, transmolprop, indexmod,indexslo,pentes_simu,modules_simu
killwaves/Z uniktest,priosep,priopoint,priolevel,priolevel,prioindex,prioslo,priointercept, rapportslo, interceptrans
killwaves/Z xxtrans, checkslo, slopematch,slopes,modules,intercepts,diffdm,diffm, diffint, cropmass,cropdefmass
killwaves/Z criblepente, slopemol, unikmod, uniksep, uniklvl,priomod, transient0, wavesourse,size0, FTMS, costfunction
killwaves/Z parent, indexintensite,intensiteold, masseold,adja, calcerror
end

function crystalinit()
wave/Z slopetrix, modules, slopes, molecule
wave/T slopenames
svar nommol
variable k,i,n,m
n=dimsize(slopetrix,0)
make/O/B/N=(n,n) cristcheck
duplicate/O molecule moletest
cristcheck=1
for(k=0;k<1;k+=1)
	for(i=k+1;i<n;i+=1)
		attribcombi(modules[i][k],slopes[i][k])
		slopetrix[i][k][]=molecule[z]
		slopenames[i][k]=nommol
		doupdate
	endfor
endfor
end

function crystalpropa()// pas bon
wave/Z slopetrix, modules, slopes, molecule
wave/T slopenames
svar nommol
variable k,i,n,m
n=dimsize(slopetrix,0)
for(k=1;k<n;k+=1)
	for(i=k+1;i<n;i+=1)
		molecule=slopetrix[i][k-1][x]-slopetrix[k][k-1][x]
		slopetrix[i][k][]=molecule[z]
		genestringpente()
		slopenames[i][k]=nommol
		doupdate
	endfor
endfor
end

function bourrin()//calcul l'atrib relative pour chaque couple
wave/Z slopetrix, modules, slopes, molecule, calcerror,modules_simu
wave/T slopenames
svar nommol
variable k,i,n,m
n=dimsize(slopetrix,0)
make/O/B/N=(n,n) cristcheck
duplicate/O molecule moletest
cristcheck=1
for(k=0;k<n;k+=1)
	for(i=k+1;i<n;i+=1)
		attribcombi(modules[i][k],slopes[i][k])
		slopetrix[i][k][]=molecule[z]
		slopenames[i][k]=nommol
		calcerror[i][k]=modules_simu[0]
		doupdate
	endfor
endfor
end

function alternattrib()//attribue par essai erreur en se rapprochant de la bone valeur
nvar inco
string lamol
wave/Z perm, molecule,maxperm,molperm_0
variable mass1,mass2, test, k,n
duplicate/O perm rand
duplicate/O molecule lafor1 lafor2 trans
lafor1=maxperm[0]*molperm_0
n=numpnts(rand)
do
rand=(round(enoise(10)))
lafor2=lafor1
	for(k=0;k<n;k+=1)
		lamol="molperm_"+num2str(k)
		duplicate/O $lamol trans
		lafor2+=trans*rand[k]
	endfor
	duplicate/O lafor2 molecule
	mass2=genestringpente()
	if(abs(inco-mass2)<abs(inco-mass1))
		lafor1=lafor2
		duplicate/O lafor1 molecule
		mass1=genestringpente()
		doupdate
		print mass1
	endif
while(mass1!=inco)

killwaves rand trans
return mass1
end

//=========ARIZONA DREAM STOPS=================


//=====MAYLIS sauvegarde des listes d'ions STARTS===========

function singlemol2matmolbag(nom)
string nom
string nomdumolbag
nomdumolbag="Molbag_"+nom
if(itemsinlist(wavelist("Molecule_*",";","")))
concatenate/O wavelist("Molecule_*",";",""), $nomdumolbag
else
Make/O/N=0 $nomdumolbag
endif
end

function chopelacharge(nom)
string nom
wave/T list_formules
variable k,n=numpnts(list_formules),poz
string transformule, nomdescharges
nomdescharges="Chargebag_"+nom
make/O/N=(n) lescharges
	for (k=0;k<n;k+=1)
		transformule=list_formules[k]
		poz=strsearch(transformule,"\S",0)
		if(poz>-1)
			lescharges[k]=str2num(transformule[poz+3]+transformule[poz+2])
		else
			lescharges[k]=0
		endif
	endfor
duplicate/O lescharges, $nomdescharges
killwaves lescharges
end

function sauvelaliste(nom)
string nom
wave/T list_formules
string nomdelaliste
nomdelaliste="Listbag_"+nom
duplicate/O list_formules $nomdelaliste
end

//ajout du contenu d'un molbag dans la liste courante

function addthemolbag(nom)
string nom
string lemolbag, lacharge, lacible
variable t=ticks
lemolbag="Molbag_"+nom
lacharge="Chargebag_"+nom
lacible="Ciblebag_"+nom
wave/Z molecule, list_targets, localcible=$lacible,simu_proba, incoy
nvar charge,krikriter
variable k,n, tranzcharge=charge, xeuma
duplicate/O $lemolbag transbag
duplicate/O $lacharge transcharge
duplicate/O molecule transmol
n=numpnts(transcharge)
for(k=0;k<n;k+=1)
	tuelesspec()
	molecule=transbag[x][k]
	charge=transcharge[k]
	genestringmol()
	simu_proba*=localcible[k][2]/incoy[0]
	addthissimu()
	list_targets[dimsize(list_targets,0)][]=localcible[k][y]
endfor
molecule=transmol[x]
genestringmol()
charge=tranzcharge
killwaves transbag, transcharge//, transmol
majliste()
majlegende()
end

//comparer deux molbags

function comparelistbag(nom1,nom2)
string nom1,nom2
string list1,list2,nomlistecom, nomlisteonly2, nomlisteonly1,molbag1,molbag2,nommolbagcom,nommolbagonly2,nommolbagonly1,chargebag1,chargebag2,nomchargebagcom,nomchargebagonly2,nomchargebagonly1
nomlistecom="Listbag_"+nom1+"_ET_"+nom2
nomlisteonly2="Listbag_only"+nom2
nomlisteonly1="Listbag_only"+nom1
nommolbagcom="Molbag_"+nom1+"_ET_"+nom2
nommolbagonly2="Molbag_only"+nom2
nommolbagonly1="Molbag_only"+nom1
nomchargebagcom="Chargebag_"+nom1+"_ET_"+nom2
nomchargebagonly2="Chargebag_only"+nom2
nomchargebagonly1="Chargebag_only"+nom1
list1="Listbag_"+nom1
list2="Listbag_"+nom2
molbag1="Molbag_"+nom1
molbag2="Molbag_"+nom2
chargebag1="Chargebag_"+nom1
chargebag2="Chargebag_"+nom2
duplicate/O/T $list1 trans1
duplicate/O/T $list2 trans2
duplicate/O $molbag1 transmolbag1
duplicate/O $molbag2 transmolbag2
duplicate/O $chargebag1 transchargebag1
duplicate/O $chargebag2 transchargebag2
variable n=numpnts(trans1), m=numpnts(trans2),k,p
make/O/N=(n,m) compmat
compmat=stringmatch(trans1[x],trans2[y])
Matrixop/O match2=sumcols(compmat)^t
Matrixop/O match1=sumcols(compmat^t)^t

make/O/N=114 transmol
make/O/T/N=0 Listbagcom, Listbagonly2, Listbagonly1
make/O/N=0 molbagcom,molbagonly1,molbagonly2
make/O/N=0 chargebagcom,chargebagonly1,chargebagonly2
for (k=0;k<m;k+=1)
transmol=0
	if (match2[k]==1)// création d'une liste des molécules communes aux deux "bag"
		insertpoints numpnts(Listbagcom),1,Listbagcom
		listbagcom[numpnts(Listbagcom)-1]=trans2[k]
		transmol=transmolbag2[x][k]
		concatenate {transmol},molbagcom
		insertpoints numpnts(chargebagcom),1,chargebagcom
		chargebagcom[numpnts(chargebagcom)-1]=transchargebag2[k]
	else // création d'une liste des molécules que possède uniquement la liste 2
		insertpoints numpnts(listbagonly2),1,Listbagonly2
		listbagonly2[numpnts(listbagonly2)-1]=trans2[k]
		transmol=transmolbag2[x][k]
		concatenate {transmol},molbagonly2
		insertpoints numpnts(chargebagonly2),1,chargebagonly2
		chargebagonly2[numpnts(chargebagonly2)-1]=transchargebag2[k]
	endif
endfor

for (k=0;k<n;k+=1)
	if (match1[k]!=1) // création d'une liste des molécules que possède uniquement la liste 1
		insertpoints numpnts(listbagonly1),1,Listbagonly1
		listbagonly1[numpnts(listbagonly1)-1]=trans1[k]
		transmol=transmolbag1[x][k]
		concatenate {transmol},molbagonly1
		insertpoints numpnts(chargebagonly1),1,chargebagonly1
		chargebagonly1[numpnts(chargebagonly1)-1]=transchargebag1[k]
	endif
endfor

duplicate/O Listbagcom $nomlistecom
duplicate/O Listbagonly2 $nomlisteonly2
duplicate/O Listbagonly1 $nomlisteonly1
duplicate/O molbagcom $nommolbagcom
duplicate/O molbagonly2 $nommolbagonly2
duplicate/O molbagonly1 $nommolbagonly1
duplicate/O chargebagcom $nomchargebagcom
duplicate/O chargebagonly2 $nomchargebagonly2
duplicate/O chargebagonly1 $nomchargebagonly1

killwaves trans1, trans2, Listbagcom, Listbagonly2, Listbagonly1,compmat,match2,match1,transmol,molbagcom,transmolbag1,transmolbag2,molbagonly1,molbagonly2,chargebagcom,transchargebag1,transchargebag2,chargebagonly1,chargebagonly2
end

//======MAYLIS sauvergarde des listes d'ion ENDS=============




Window Graph0() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(594.75,369.5,823.5,571.25) histoiso
	AppendToGraph/T histoiso
	AppendToGraph/B=top2 histoiso
	AppendToGraph/B=top3 histoiso
	ModifyGraph userticks(top)={defmassreliso12,elements}
	ModifyGraph userticks(top2)={defmassreliso13,elements}
	ModifyGraph userticks(top3)={defmassreliso23,elements}
	ModifyGraph mode=5
	ModifyGraph rgb=(0,0,0)
	ModifyGraph grid(top)=1,grid(top2)=1,grid(top3)=1
	ModifyGraph tick(left)=3
	ModifyGraph noLabel(left)=2
	ModifyGraph lblMargin(top)=32
	ModifyGraph axOffset(left)=-6,axOffset(bottom)=-1,axOffset(top)=0.5
	ModifyGraph axThick(left)=0
	ModifyGraph gridRGB(top2)=(0,52224,0),gridRGB(top3)=(65280,0,0)
	ModifyGraph axRGB(top)=(24576,24576,65280),axRGB(top2)=(0,52224,0),axRGB(top3)=(65280,0,0)
	ModifyGraph tlblRGB(top)=(24576,24576,65280),tlblRGB(top2)=(0,52224,0),tlblRGB(top3)=(65280,0,0)
	ModifyGraph alblRGB(top)=(24576,24576,65280),alblRGB(top2)=(0,52224,0),alblRGB(top3)=(65280,0,0)
	ModifyGraph fStyle(top)=1, fStyle(top2)=1, fStyle(top3)=1
	ModifyGraph lblPos(bottom)=25
	ModifyGraph freePos(top2)=-109
	ModifyGraph freePos(top3)=-81
	Label top "\\Z061\\Sst\\M\\Z06 - 2\\Snd"
	Label top2 "\\Z061\\Sst\\M\\Z06 - 3\\Srd"
	Label top3 "\\Z062\\Snd\\M\\Z06 - 3\\Srd"
	Button button0,pos={185,3},size={50,20},proc=ButtonZ2,title="Zoom2"
	Button button1,pos={57,3},size={50,20},proc=ButtonZ1,title="Zoom1"
	Button button2,pos={121,3},size={50,20},proc=ButtonZ0,title="Reset"
EndMacro

//==========MANUALLY ENTERED MASS MODIFICATION STARTS======
// 12 april 2012
// allows to calculate on manually entered mass in the current m/z field
// Find in data check box is created in attributeur and has action on agregateur, dmvm and attributeur

Function SetVarProcInco(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
nvar inco,findindata
wave/Z wave1stripped,wave0stripped, incox, incoxdm,incoy
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
				if(findindata==1)
					incox=wave0stripped(prochepic(inco,0))
					incoy=wave1stripped(prochepic(inco,0))
					inco=wave0stripped(prochepic(inco,0))
					setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
					setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
				else
					incox=inco
				endif
			incoxdm=defmass(incox)
			majMarkerFindOrNot()
			break
	endswitch
	return 0
End

Function CheckFindInData_proc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
nvar inco
wave/Z wave1stripped,wave0stripped, incox, incoxdm,incoy
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			majMarkerFindOrNot()
			if(checked==1)
					incox=wave0stripped(prochepic(inco,0))
					incoy=wave1stripped(prochepic(inco,0))
					inco=wave0stripped(prochepic(inco,0))
					setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
					setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
				else
					incox=inco
				endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function majMarkerFindOrNot()
nvar findindata,inco
wave/Z wave1stripped,incoy
variable k=wavemax(wave1stripped)/2
if(findindata==1)
	ModifyGraph/W=agregateur opaque(incoy)=1,marker(incoy)=42,msize(incoy)=5
	ModifyGraph/W=dmvm opaque(incoxdm)=1,marker(incoxdm)=42,msize(incoxdm)=5
	ErrorBars/W=agregateur incoy OFF
elseif(findindata==0)
	ModifyGraph/W=agregateur marker(incoy)=28
	ModifyGraph/W=dmvm marker(incoxdm)=28
	ErrorBars/W=agregateur/T=0 incoy Y,const=k
	incoy=k
endif
end

//==========MANUALLY ENTERED MASS MODIFICATION ENDS======

//==========Visual modes switcher STARTS====================
//12 avril 2012

function makeseg()
wave/Z incox, incoy, incoxdm, dataindex,wave1stripped, diffm,diffdm
variable k
variable/G twosig
findvalue/V=(incoy[0]) wave1stripped
findvalue/V=(v_value) dataindex
k=v_value
if(k<= twosig)
make/O/N=(3*numpnts(wave1stripped)) segzm, segzdm
segzm=(mod(x,3)==0)*incox+(mod(x-1,3)==0)*(diffm[x/3][k]+incox)
segzdm=(mod(x,3)==0)*incoxdm+(mod(x-1,3)==0)*(diffdm[x/3][k]+incoxdm)
segzm*=segzm/segzm
segzdm*=segzdm/segzdm
else
make/O/N=(0) segzm, segzdm
endif
end

//==========Visual modes switcher STARTS====================



////rentree2013

function traceresneeded(b1,b2,pas)
variable b1,b2,pas
nvar inco
nvar krikriter
wave/Z list_molperm, molecule,simu_mass
string lawave
variable k, n=(max(b1,b2)-min(b2,b1))/pas, i=1
make/o/n=(n) mass2attrib, attribmass, resneed
mass2attrib=b1
attribmass=b1
resneed=0
krikriter=1
display resneed vs mass2attrib
ModifyGraph log(left)=1
appendtograph resneed vs attribmass
ModifyGraph rgb(resneed#1)=(0,15872,65280)
ModifyGraph mode(resneed#1)=3
ModifyGraph marker(resneed#1)=9
	for(k=0;k<n;k+=1)
		inco=b1+k*pas
		mass2attrib[k]=inco
		analysemass("arg")
		setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
		setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
		lawave=stringfromlist(0,wavelist("molprop_*",";",""))
		molecule=0
		tuelesspec()
		duplicate/O $lawave molecule
		genestringmol()
		wave/Z deltappm
			do
				resneed[k]=abs(deltappm[0]-deltappm[i])
				i+=1
			while(resneed[k]==0)
		attribmass[k]=simu_mass[0]
		doupdate
	endfor
end

function spandensite(b1,b2)
variable b1,b2
nvar inco
nvar krikriter
wave/Z list_molperm, molecule,simu_mass,masseclassee,deltappm,lesmassescriblees
string lawave
variable k, pas,i=0
make/o/n=0 attribmass, resneed, taillecrible
krikriter=1
display resneed vs attribmass
appendtograph taillecrible vs attribmass
inco=b1
	do
		analysemass("arg")
		setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
		setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
		lawave=stringfromlist(0,wavelist("molprop_*",";",""))
		molecule=0
		tuelesspec()
		duplicate/O $lawave molecule
		genestringmol()
		insertpoints 0,1, attribmass, resneed, taillecrible
		attribmass[0]=masseclassee[0]
			do
				i+=1
				resneed[0]=deltappm[i]
			while(resneed[0]<=0)
		taillecrible[0]=numpnts(lesmassescriblees)
		doupdate
		inco=masseclassee[i]
	while(inco<b2)
end

function cleanmass()
wave/Z attribmass
variable k, n=numpnts(attribmass)
make/o/N=1 cleanedattribmass
print n
cleanedattribmass[0]=attribmass[0]
	for(k=1;k<n;k+=1)
		if(attribmass[k-1]!=attribmass[k])
			insertpoints 0,1, cleanedattribmass
			cleanedattribmass[0]=attribmass[k]
		endif
	endfor
sort cleanedattribmass, cleanedattribmass
print numpnts(cleanedattribmass)
end

function historework()
wave/Z slopes
make/o/n=(numpnts(slopes)) slopelin
slopelin=slopes[x]
wavetransform zapnans slopelin
sort slopelin slopelin
Differentiate slopelin/D=slopelin_diff
findlevels slopelin_diff 0
wave/Z W_FindLevels
duplicate/o W_FindLevels slopcrop
slopcrop=slopelin[W_FindLevels[x]]
variable k,n=numpnts(slopcrop)
make/o/n=1 slopunik, slopunikcount
slopunik[0]=slopcrop[0]
slopunikcount[0]=1
for(k=1;k<n;k+=1)
	if(slopcrop[k]!=slopcrop[k-1])
		insertpoints 0,1, slopunik, slopunikcount
		slopunik[0]=slopcrop[k]
		slopunikcount[0]=1
	else
		slopunikcount[0]+=1
	endif
endfor
sort slopunik, slopunikcount, slopunik
end




Function Proc_setintenscumthres(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
nvar twosig, intenscumthres
wave/Z dataindex, logwave1
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
	endswitch
	twosig=intenstat(intenscumthres)
	valse(twosig)
	variable coul=logwave1[dataindex[intenstat(intenscumthres)]]
	ModifyGraph/W=dmvm zColor(wave2stripped)={logwave1,coul ,*,YellowHot,0};
	ModifyGraph/W=dmvm zColorMin(wave2stripped)=(65535,65535,65535)
	ModifyGraph/W=agregateur zColor(wave1stripped)={logwave1,coul ,*,YellowHot,0};
	ModifyGraph/W=agregateur zColorMin(wave1stripped)=(65535,65535,65535)
	calcinfo()
	return 0
End


function ripradio()
wave/Z wave0stripped, wave2stripped
variable j, k, n=numpnts(wave0stripped)
make/o/n=0 interceptradio, locaradio
display interceptradio
for(j=0; j<(n-1); j+=1)
	for(k=j+1;k<n;k+=1)
		if(	(	(wave2stripped[k]-wave2stripped[j])/(wave0stripped[k]-wave0stripped[j])	)==1 	)
			insertpoints 0,1, interceptradio, locaradio
			interceptradio[0]=wave2stripped[j]- ((wave2stripped[k]-wave2stripped[j])/(wave0stripped[k]-wave0stripped[j])) * wave0stripped[j]
			locaradio[0]=j
		endif
	endfor
	doupdate
endfor
end

function ripradio2()
wave/Z wave0stripped, wave2stripped, wave1stripped
variable j, k, n=numpnts(wave0stripped), teupen, intercept
make/o/n=0 interceptradio, locaradio, countradio
display interceptradio
for(j=0; j<(n-1); j+=1)
	k=j+1
	teupen=(wave2stripped[k]-wave2stripped[j])/(wave0stripped[k]-wave0stripped[j])
		if(teupen==1)
			intercept=wave2stripped[j]- teupen * wave0stripped[j]
			insertpoints 0,1, interceptradio, locaradio, countradio
			interceptradio[0]=wave2stripped[j]- ((wave2stripped[k]-wave2stripped[j])/(wave0stripped[k]-wave0stripped[j])) * wave0stripped[j]
			locaradio[0]=j
			countradio[0]=1/abs(wave1stripped[j]-wave1stripped[k])
		endif
	doupdate
	make/o histo
	histogram/B={wavemax(interctradio),-1,-wavemin(interceptradio)+wavemax(interctradio)} interceptradio histo
endfor
end

function radio2roi(seuil)
variable seuil
variable moy
wave/Z wave0stripped, wave1stripped, wave2stripped
duplicate/o wave0stripped test, roipnts
test=round(wave0stripped)
make/o histo
histogram/B={wavemin(test),1, wavemax(test)-wavemin(test)} test histo
duplicate/o histo, histo_int, fit_histo
integrate histo_int
integrate histo_int
CurveFit/NTHR=0/Q poly_XOffset 10,  histo_int /D=fit_histo
differentiate fit_histo
differentiate fit_histo
moy=mean(histo)
histo-=fit_histo+seuil*moy
roipnts=(histo(test)>0)*x
roipnts=roipnts/roipnts+roipnts-1
wavetransform zapnans roipnts
duplicate/o roipnts roi1,roi0,roi2
roi1=wave1stripped[roipnts[x]] ; roi0=wave0stripped[roipnts[x]] ; roi2=wave2stripped[roipnts[x]]
killwaves/Z test
end

function roi2data(newname)
string newname
string les1, les0
les1=newname+"_1stripped"
les0=newname+"_0stripped"
duplicate/O roi1 $les1
duplicate/O roi0 $les0
if(waveexists(roi3))
	string les3=newname+"_3stripped"
	duplicate/O roi3 $les3
endif
listelesdonnees()
end

function invertroi()
wave/Z roipnts, wave1stripped, wave2stripped, wave0stripped, roi1, roi2, roi0
sort/r roipnts roipnts
duplicate/o roipnts oldroipnts
variable k, n=numpnts(oldroipnts)
make/o/n=(numpnts(wave0stripped)) roipnts
roipnts=x
for(k=0;k<n;k+=1)
	deletepoints oldroipnts[k],1,roipnts
endfor
killwaves oldroipnts
roipnts2roiz()
calcinfo()
end

function roipnts2roiz()
wave/Z roipnts, wave1stripped, wave2stripped, wave0stripped
duplicate/o roipnts, roi1, roi2, roi0
roi0=wave0stripped[roipnts[x]]
roi1=wave1stripped[roipnts[x]]
roi2=wave2stripped[roipnts[x]]
sort roi0, roi0, roi1,roi2
if(waveexists(wave3stripped))
	wave wave3stripped
	duplicate/O roipnts, roi3
	roi3=wave3stripped[roipnts[x]]
	sort roi0, roi3
endif
end

function calcinfo()
variable t=ticks
nvar logscaleon, intenscumthres, livedatastat
if(livedatastat)
make/o histopoints, histoint, historoi
variable b1,b2, xeu,b3,b4,coul,tot
wave/Z wave1stripped, logwave1, cursintx, dataindex, roi1
string wavetoanal="wave1stripped", roi2anal="roi1"
getaxis/q/W=agregateur left
if(logscaleon==1)
	wavetoanal="logwave1"
	duplicate/o $roi2anal logroi1
	duplicate/o wave1stripped logwave1
	multithread logwave1=log(wave1stripped)
	roi2anal="logroi1"
	logroi1=log(roi1)
	b1=log(v_min)
	b2=log(v_max)
elseif(logscaleon==0)
	b1=v_min
	b2=v_max
endif
getaxis/q/W=agregateur bottom
b3=(mass2pnt(v_min))
b4=(mass2pnt(v_max))
histogram/B={b1, (b2-b1)/100, 101}/R=[b3,b4] $wavetoanal histopoints
histogram/B={b1, (b2-b1)/100, 101}/R=[b3,b4]/W=wave1stripped $wavetoanal histoint
if(numpnts(roi1)!=0)
	histogram/B={b1, (b2-b1)/100, 101}/W=roi1 $roi2anal historoi
elseif(numpnts(roi1)==0)
	historoi=0
endif
duplicate/o histoint histocumint, histoabcisse, histocumbar, histocumpoints, histocolorint
multithread histoabcisse=x
multithread histocolorint=x
if(logscaleon==1)
	multithread histoabcisse=10^(x)
	coul=logwave1[dataindex[intenstat(intenscumthres)]]
	ModifyGraph/W=datastat zColor(histoabcisse#2)={histocolorint,coul,*,YellowHot,0}
	modifygraph/w=datastat log(left)=1
	setaxis/W=datastat left, 10^b1, 10^b2
else
	coul=wave1stripped[dataindex[intenstat(intenscumthres)]]
	ModifyGraph/W=datastat zColor(histoabcisse#2)={histocolorint,coul,*,YellowHot,0}
	modifygraph/w=datastat log(left)=0
	setaxis/W=datastat left, b1,b2
endif
	coul=logwave1[dataindex[intenstat(intenscumthres)]]
	ModifyGraph/W=dmvm zColor(wave2stripped)={logwave1,coul ,*,YellowHot,0};
	ModifyGraph/W=dmvm zColorMin(wave2stripped)=(65535,65535,65535)
	ModifyGraph/W=agregateur zColor(wave1stripped)={logwave1,coul ,*,YellowHot,0};
	ModifyGraph/W=agregateur zColorMin(wave1stripped)=(65535,65535,65535)
multithread histocumpoints=sum(histopoints,x,b2)
multithread histocumint=sum(histoint,x,b2)
xeu=1/wavemax(histoint)
tot=1/sum(wave1stripped)
fastop histocumint=(tot)*histocumint
fastop histoint=(xeu)*histoint
fastop historoi=(xeu)*historoi
fastop histocumbar=1-histocumint
//print b1, b2, (b2-b1)/100
cursintx[0]=wavemax(histocumpoints)/2
killwaves/Z logroi1
endif
//print "calcinfo :", ticks-t //1 ou 2 ticks
end

threadsafe function mass2pnt(m)
variable m
variable p
wave/Z wave0stripped
findlevel/q wave0stripped m
p= v_levelx
if(m>wavemax(wave0stripped))
	p=numpnts(wave0stripped)-1
elseif(m<wavemin(wave0stripped))
	p=0
endif
return p
end

Function CheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
nvar logscaleON

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			if(logscaleON==1)
			ModifyGraph/W=agregateur log(left)=1
			ModifyGraph/W=elaborateur log(left)=1
			else
			ModifyGraph/W=agregateur log(left)=0
			ModifyGraph/W=elaborateur log(left)=0
			endif
			break
	endswitch
calcinfo()
	return 0
End

Window datastat() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(830.25,42.5,939.75,484.25)/T simu_proba vs simu_plot as "datastat"
	AppendToGraph histoabcisse vs histopoints
	AppendToGraph histoabcisse vs histocumpoints
	AppendToGraph/T histoabcisse vs histocumint
	AppendToGraph/T histoabcisse vs histoint
	AppendToGraph cursint vs cursintx
	AppendToGraph/T cursint vs cumintx
	AppendToGraph cursint vs cumpointsx
	AppendToGraph/T histoabcisse vs historoi
	AppendToGraph incoy[0,0] vs cursintx[0,0]
	AppendToGraph/T simu_proba vs simu_plot
	ModifyGraph userticks(top)={cumintx,labcumintx}
	ModifyGraph userticks(bottom)={cumpointsx,labcumpointsx}
	ModifyGraph gbRGB=(47872,47872,47872)
	ModifyGraph mode(simu_proba)=3,mode(cursint)=4,mode(cursint#1)=3,mode(cursint#2)=3
	ModifyGraph mode(incoy)=3,mode(simu_proba#1)=3
	ModifyGraph marker(simu_proba)=18,marker(histoabcisse)=10,marker(histoabcisse#3)=10
	ModifyGraph marker(cursint)=25,marker(cursint#1)=41,marker(cursint#2)=42,marker(incoy)=42
	ModifyGraph lSize(histoabcisse#4)=2
	ModifyGraph rgb(simu_proba)=(0,52224,52224),rgb(histoabcisse#2)=(0,0,0),rgb(histoabcisse#3)=(0,0,65280)
	ModifyGraph rgb(cursint)=(26368,0,52224),rgb(cursint#1)=(0,0,0),rgb(histoabcisse#4)=(26368,0,52224)
	ModifyGraph rgb(incoy)=(26368,0,52224),rgb(simu_proba#1)=(0,0,0)
	ModifyGraph msize(cursint)=5,msize(simu_proba#1)=2
	ModifyGraph opaque(cursint)=1,opaque(incoy)=1
	ModifyGraph useMrkStrokeRGB(simu_proba)=1
	ModifyGraph mrkStrokeRGB(simu_proba)=(65535,65535,65535)
	ModifyGraph zColor(histoabcisse#2)={histocolorint,6.73638861318596,*,YellowHot}
	ModifyGraph zColorMin(histoabcisse#2)=(65535,65535,65535)
	ModifyGraph textMarker(simu_proba#1)={simu_proba,"default",0,0,6,5.00,0.00}
	ModifyGraph log(left)=1
	ModifyGraph fStyle(left)=1
	ModifyGraph lblMargin(left)=3,lblMargin(top)=1
	ModifyGraph axOffset(left)=-3.85714,axOffset(top)=-0.266667,axOffset(bottom)=-0.266667
	ModifyGraph axThick(left)=2
	ModifyGraph lblLatPos(left)=-105,lblLatPos(top)=-7
	Label left "Intensity"
	Label top "Intensity CDF"
	Label bottom "Points CDF"
	SetAxis left 20518.76,5490769.5
	ErrorBars/T=0 histoabcisse X,wave=(,histopoints)
	ErrorBars/T=0 histoabcisse#1 X,wave=(,histocumpoints)
	ErrorBars/T=0 histoabcisse#2 X,wave=(histocumbar,)
	ErrorBars/T=0 histoabcisse#3 X,wave=(,histoint)
	ErrorBars cursint X,wave=(cursintx,cursintx)
	ErrorBars/T=0 cursint#1 Y,wave=(ebarcmpntsx,)
	ErrorBars/T=0 cursint#2 Y,wave=(,cursint)
	ControlBar 40
	SetVariable setvar0,pos={0,287},size={50,16},proc=SetVarProc_4,title=" "
	SetVariable setvar0,limits={-inf,inf,0},value= cursint[0]
	SetVariable setvar1,pos={0,287},size={50,16},disable=1,proc=SetVarProc_5,title=" "
	SetVariable setvar1,limits={-inf,inf,0},value= cursint[1]
	CheckBox check0,pos={3,3},size={74,14},title="Live update"
	CheckBox check0,variable= livedatastat
	Button callHough,pos={2,18},size={89,20},proc=HoughTcaller,title="HoughT Analysis"
	Button callHough,fColor=(65535,65535,65535)
	SetWindow kwTopWin,hook=getdatastat,hookevents=7
EndMacro


function majcursint(val)
variable val
nvar logscaleon
wave/Z cursint, histocumint, histocumpoints, histoabcisse,cursintx
cursint[0]=val
duplicate/o cursint, cumintx, cumpointsx, ebarcmpntsx
if(logscaleon)
	cumintx=histocumint(log(cursint[x]))
	cumpointsx=round(histocumpoints(log(cursint[x])))
else
	cumintx=histocumint((cursint[x]))
	cumpointsx=round(histocumpoints((cursint[x])))
endif
ebarcmpntsx=wavemax(histoabcisse)-cursint
make/t/o/n=(numpnts(cursint)) labcumintx, labcumpointsx
labcumintx=num2str(cumintx)
labcumpointsx=num2str(cumpointsx)
cursintx=wavemax(cursintx)
setvariable setvar0 win=datastat, title=" ", pos={0,pixelfromaxisval("datastat","left",cursint[0])}, value=cursint[0], proc=SetVarProc_4, frame=1, noedit=0, limits={-inf,inf,0}
setvariable setvar1 win=datastat, title=" ", pos={0,pixelfromaxisval("datastat","left",cursint[1])}, value=cursint[1], proc=SetVarProc_5, frame=1, noedit=0, limits={-inf,inf,0}
	if(cursint[0]==cursint[1])//numpnts(cursint)==2 && 
		setvariable setvar1,win=datastat, disable=1
	else
		setvariable setvar1,win=datastat, disable=0
	endif
end

function getdatastat(info)
string info
make/o/d/n=1 pozint
pozint=axisvalfrompixel("datastat","left",str2num(stringbykey("mousey",info)))
wave/Z cumpointsx,dataindex
if(stringmatch(stringbykey("event",info),"mousemoved")==1 && str2num(stringbykey("modifiers",info))==1)
majcursint(pozint[0])
roifromcursint()
endif
if(stringmatch(stringbykey("event",info),"mousedown")==1)
make/o/d/n=1 cursint,cursintx
majcursint(pozint[0])	
insertpoints 0,1, cursint, cursintx
majcursint(pozint[0])	
endif
if(stringmatch(stringbykey("event",info),"mouseup")==1)
majcursint(pozint[0])	
	if(numpnts(cursint)==2 && cursint[0]==cursint[1])
		deletepoints 0,1, cursint, cursintx
	endif
	roifromcursint()
endif
return 1
end

function roifromcursint()
wave/Z cursint,dataindex, sortedint,wave0stripped,wave1stripped
//make/o/n=((abs(cumpointsx[1]-cumpointsx[0])+1)*(cumpointsx[1]-cumpointsx[0]!=0)) roipnts
//roipnts=dataindex[x+min(cumpointsx[1],cumpointsx[0])-1]
//roipnts2roiz()
cropbetween(cursint[0],cursint[1],wave0stripped,wave1stripped)
duplicate/O outy roi1, roi2
duplicate/O outx roi0, roipnts
killwaves outy, outx
roi2=defmass(roi0)
sort roi0, roi0, roi1, roi2
roipnts=mass2pnt(roi0)
calcinfo()
end

function seekheavy()
string leswaves, lawave
leswaves=wavelist("*",";","")
variable k, n= itemsinlist(leswaves)
make/o/t/n=(n) wavesnames
make/o/d/n=(n) wavesnbpnt
for(k=0;k<n;k+=1)
	lawave=stringfromlist(k,leswaves)
	wavesnames[k]=lawave
	wavesnbpnt[k]=numpnts($lawave)
endfor
sort/R wavesnbpnt, wavesnbpnt, wavesnames
edit wavesnames, wavesnbpnt
end

function exaustoslo(nbpnt)
variable nbpnt
wave/Z wave2stripped, wave1stripped, wave0stripped, histoslo
variable k, j, n=numpnts(wave2stripped), laslo
make/o exaustivehistoslo
histogram/B={-0.01,0.02/nbpnt,nbpnt} wave1stripped exaustivehistoslo
duplicate/o exaustivehistoslo histoslo
exaustivehistoslo=0
for(k=1;k<n;k+=1)	
	for(j=0;j<k;j+=1)
		laslo=(wave2stripped[j]-wave2stripped[k])/(wave0stripped[j]-wave0stripped[k])
		if(abs(laslo)<=0.01)
 		exaustivehistoslo[x2pnt(exaustivehistoslo,laslo)]+=1
 		endif
	endfor
endfor
histoslo=exaustivehistoslo
histo2rose()
doupdate
end

function histo2rose()
wave/Z histoslo, wave0stripped
variable maxx, maxy, contrainte
duplicate/O histoslo xplode yplode
xplode=(histoslo)*cos(atan(pnt2x(histoslo,p)))
yplode=(histoslo)*sin(atan(pnt2x(histoslo,p)))
maxx=max(wavemax(xplode), abs(wavemin(xplode)))
maxy=max(wavemax(yplode), abs(wavemin(yplode)))
contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
xplode*=contrainte
yplode*=contrainte
xplode+=wavemin(wave0stripped)
end

////RELEASE post-rentree2013

//---------BEGIN MAJ Mai 2013

function updateDistanceDMVM(xx,yy)
variable xx,yy
nvar inco
wave/Z Wave0stripped, wave2stripped, incox,incoy,incoxdm, wave1stripped
variable stretchx, stretchy
duplicate/O wave0stripped distance
getaxis/Q/W=dmvm bottom
stretchx=v_max-v_min
getaxis/Q/W=dmvm left
stretchy=v_max-v_min
distance=( (wave0stripped-xx)/stretchx )^2 + ( (wave2stripped-yy)/stretchy )^2
findlevel/Q distance (wavemin(distance))
return v_levelx
end

///----modif getdmvmMSMSmode
///----modif attribcombi
///---------END MAJ Mai 2013


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////
///// Fusion percée 1.07_dev de Roland et 1.08_dev de FROD
/////-fullauto updated
/////-partie Djevahirdjian Leo 2013 annule et remplace Djevahirdjian mai 2012
/////
///////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////  REPRESENTATION TOOLS  July 2013 Roland ///////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function extractreport ()
wave/T titres
wave/T rapport
variable numtitre = dimsize(titres,1), k
variable taille=dimsize (rapport,0)
string nomwav	
make/O/N=(taille-1) element_O, element_H, element_N, element_C
element_O=0
element_H=0
element_C=0
element_N=0
for (k=0; k<numtitre; k+=1)
	if (k<4)
		nomwav=titres[0][k]
	else
		nomwav="element_"+titres[0][k]
	endif
	make/O/N=(taille) $nomwav
	wave w= $nomwav
	w=str2num(Rapport[x][k+1])
	deletepoints 0,1, w
endfor
end

Function appendelements(list)
	String list					
	String theWave
	Variable index=0

	do
		theWave = StringFromList(index, list)
		if (strlen(theWave) == 0)
			break							
		endif
		AppendToGraph $theWave vs MassInconnue
		index += 1
	while (1)							
End


Window compare_attributions() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(490.5,37.25,1433.25,777.5)/L=qual Proximite vs MassInconnue
	DoWindow/T compare_attributions,"compare_attributions"
	Label qual "proximité"
	AppendToGraph/L=delt Delta vs MassInconnue
	Label delt "Delta ppm"
	AppendToGraph element_H vs MassInconnue
	AppendToGraph element_C vs MassInconnue
	AppendToGraph element_N vs MassInconnue
	AppendToGraph element_O vs MassInconnue
	Label left "Nombre d'atomes"
	ModifyGraph mode(Proximite)=1,mode(Delta)=1,mode(element_H)=3,mode(element_C)=3
	ModifyGraph mode(element_N)=3,mode(element_O)=3
	ModifyGraph marker(element_H)=41,marker(element_C)=41,marker(element_N)=41,marker(element_O)=41
	ModifyGraph lSize(element_H)=3,lSize(element_C)=3,lSize(element_N)=3,lSize(element_O)=3
	ModifyGraph rgb(Proximite)=(0,0,0),rgb(Delta)=(0,0,65280),rgb(element_H)=(65280,43520,0)
	ModifyGraph rgb(element_C)=(0,0,0),rgb(element_N)=(0,0,65280),rgb(element_O)=(65280,0,0)
	ModifyGraph msize(element_H)=3,msize(element_C)=3,msize(element_N)=3,msize(element_O)=3
	ModifyGraph mrkThick(element_H)=1.5,mrkThick(element_C)=1.5,mrkThick(element_N)=1.5
	ModifyGraph mrkThick(element_O)=1.5
	ModifyGraph zmrkSize(element_H)={Intensitepic,*,*,1,10}
	ModifyGraph grid(delt)=1
	ModifyGraph log(qual)=1
	ModifyGraph axOffset(left)=-0.222222
	ModifyGraph lblPos(left)=51
	ModifyGraph freePos(qual)=0
	ModifyGraph freePos(delt)=0
	ModifyGraph axisEnab(qual)={0.6,0.8}
	ModifyGraph axisEnab(delt)={0.8,1}
	ModifyGraph axisEnab(left)={0,0.6}
	ModifyGraph lblPos(qual)=84,lblPos(delt)=84,lblPos(left)=78
	Cursor/P A Proximite 700;Cursor/P B Delta 700;Cursor/P/S=2 C element_H 700;Cursor/P/S=2 D element_C 700;Cursor/P/S=2 E element_N 700;Cursor/P/S=2 F element_O 700
	Legend/C/N=text0/J/F=0/A=MC/X=-42.50/Y=-0.35 "\\s(delta) deltaPPM\r\\s(proximite) proximite\r\\s(element_H)H\r\\s(element_C)C\r\\s(element_N)N\r\\s(element_O)O"
	ControlBar 31
	ModifyGraph cbRGB=(0,65280,65280)
	ValDisplay valdisp1,pos={0,2},size={100,14},title="proximite"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000},value= #"qual"
	ValDisplay valdisp2,pos={111,2},size={90,14},title="masse"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000},value= #"mas"
	ValDisplay valdisp3,pos={211,2},size={90,14},title="delta"
	ValDisplay valdisp3,limits={0,0,0},barmisc={0,1000},value= #"del"
	ValDisplay valdisp4,pos={311,2},size={80,14},title="num C"
	ValDisplay valdisp4,limits={0,0,0},barmisc={0,1000},value= #"c"
	ValDisplay valdisp5,pos={401,2},size={80,14},title="num H"
	ValDisplay valdisp5,limits={0,0,0},barmisc={0,1000},value= #"h"
	ValDisplay valdisp6,pos={491,2},size={80,14},title="num N"
	ValDisplay valdisp6,limits={0,0,0},barmisc={0,1000},value= #"n"
	ValDisplay valdisp7,pos={581,2},size={80,14},title="num O"
	ValDisplay valdisp7,limits={0,0,0},barmisc={0,1000},value= #"o"
	ValDisplay valdisp8,pos={671,2},size={110,14},title="intensite"
	ValDisplay valdisp8,limits={0,0,0},barmisc={0,1000},value= #"intpic"
	Button button1,pos={791,0},size={90,30},proc=removecursorpoint,title="kill this point"
	Button button2,pos={891,0},size={80,30},proc=dorepresentations,title="represent"
	Button button11,pos={981,00},size={35,30},proc=goprev,title="prev"
	Button button12,pos={1026,00},size={35,30},proc=gonext,title="next"
	ValDisplay valdisp9,pos={1071,1},size={90,12},title="numpics",barmisc={0,1000},value= #"numpic"
	ValDisplay valdisp10,pos={1071,16},size={90,12},title="numpnts",barmisc={0,1000},value= #"numpoint"
	Button button13,pos={1171,00},size={60,30},proc=stopTheHook,title="Hook OFF"

	Button button3,pos={0,18},size={100,12},proc=sortprox,title="sort prox"
	Button button4,pos={111,18},size={90,12},proc=sortmas,title="sort mass"
	Button button5,pos={211,18},size={90,12},proc=sortdelta,title="sor delta"
	Button button6,pos={311,18},size={80,12},proc=sortC,title="sort C"
	Button button7,pos={401,18},size={80,12},proc=sortH,title="sort H"
	Button button8,pos={491,18},size={80,12},proc=sortN,title="sort N"
	Button button9,pos={581,18},size={80,12},proc=sortO,title="sort O"
	Button button10,pos={671,18},size={110,12},proc=sortInt,title="sort int"

	SetWindow kwTopWin,hook=get_compare,hookevents=1
EndMacro

function stopTheHook(ctrlName) : ButtonControl
	String ctrlName
	SetWindow kwTopWin, hook=$""
	Button button13,proc=startTheHook,title="Hook ON"
end

function startTheHook(ctrlName) : ButtonControl
	String ctrlName
	SetWindow kwTopWin,hook=get_compare,hookevents=1
	Button button13,proc=stopTheHook,title="Hook OFF"
end

function goprev(ctrlName) : ButtonControl
	String ctrlName
	variable CsrPoint = (pcsr(A)-1)
	cursor A proximite CsrPoint
	cursor B Delta CsrPoint
	PutCsrElements(CsrPoint)
	getcursorpoint("ba")
end
function gonext(ctrlName) : ButtonControl
	String ctrlName
	variable CsrPoint = (pcsr(A)+1)
	cursor A proximite CsrPoint
	cursor B Delta CsrPoint
	PutCsrElements(CsrPoint)
	getcursorpoint("ba")
end
function sortprox(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort proximite, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortmas(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort massinconnue, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortdelta(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort delta, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortC(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort element_C, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortH(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort element_H, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortN(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort element_N, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortO(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort element_O, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end
function sortInt(ctrlName) : ButtonControl
	String ctrlName
	string elements = wavelist ("element_*", ",", "")
	execute "sort Intensitepic, "+ elements+" massinconnue, delta, Intensitepic, proximite"
end


function get_compare(info)
string info
nvar Findindata,inco
variable lepoint
make/o/d/n=1 pozx, pozy, controlenfonce
pozx=str2num(stringbykey("mousex",info))
pozy=str2num(stringbykey("mousey",info))
controlenfonce=str2num(stringbykey("MODIFIERS",info))
if(stringmatch(stringbykey("event",info),"mouseup")==1)
	if(Findindata==1)
		lepoint=str2num(stringbykey("HITPOINT",TraceFromPixel (pozx[0],pozy[0],"")))
		cursor A proximite lepoint
		cursor B Delta lepoint
		PutCsrElements(lepoint)
		if (controlenfonce(0)==8)
			removecursorpoint("bob")
		endif
		getcursorpoint("ba")
	endif
endif
return 1
end



Function PutCsrElements(v_value)
	Variable v_value
	String list=wavelist ("element_*", ";", "")					
	String theWave
	Variable index=0
	string nomcurs= "CDEFGHIJ"
	string ordre
	do
		theWave = StringFromList(index, list)
		if (strlen(theWave) == 0)
			break							
		endif
		ordre ="cursor/S=2 "+nomcurs[index] +", " + theWave +","+ num2str(v_value)
		execute ordre
		index += 1
	while (1)							
End

function getcursorpoint(ctrlName) : ButtonControl
	String ctrlName
	variable pt=pcsr(A)
	variable/G qual, mas, del, intpic, c, h, n, o,numpic, numpoint
	wave w1=proximite
	qual=w1[pt]
	ValDisplay valdisp1,value=#"qual"
	wave w2=massinconnue
	mas=w2[pt]
	ValDisplay valdisp2,value=#"mas"
	wave w3=delta
	del=w3[pt]
	ValDisplay valdisp3,value=#"del"
	wave w4=element_C
	C=w4[pt]
	ValDisplay valdisp4,value=#"c"
	wave w5=element_H
	H=w5[pt]
	ValDisplay valdisp5,value=#"h"
	wave w6=element_N
	N=w6[pt]
	ValDisplay valdisp6,value=#"n"
	wave w7=element_O
	O=w7[pt]
	ValDisplay valdisp7,value=#"o"
	wave w8=Intensitepic
	intpic=w8[pt]
	ValDisplay valdisp8,value=#"intpic"
	wave/T w9=settings
	numpic=str2num(w9[0][1])
	ValDisplay valdisp9,value=#"numpic"
	wave w10=Intensitepic
	numpoint=numpnts(w10)
	ValDisplay valdisp10,value=#"numpoint"
	
end



function dorepresentations(ctrlName) : ButtonControl
	String ctrlName
	stopTheHook("bob")
	variable pt=pcsr(A)
	wave w1=proximite
	wave w2=massinconnue
	wave w3=intensitepic
	wave w4=delta
	wave w5=element_H
	wave w6=element_C
	wave w7=element_N
	wave w8=element_O
	sort  w8, w1, w2, w3, w4, w5, w6, w7, w8
	Findlevel /Q w8, 1
	if (V_flag==0)
		Make/O/N=(V_levelX) hsurc, nsurc, massShort, methyleneexcess,carbonexcess,composite, DBE
	else
		duplicate/O w8  hsurc, nsurc, massShort, methyleneexcess,carbonexcess,composite, DBE
	endif
	duplicate/O w8 zubarev
	duplicate/O w8 OsurC, HsurClong
	massShort=w2(x)
	hsurc=((w5(x)-1)-w7(x))/w6(x)
	HsurClong=((w5-1)-w7)/w6
	NsurC=w7(x)/w6(x)
	OsurC=w8/w6
	DBE=w6(x)-(w5(x)-1)/2+w7(x)/2+1
	methyleneexcess=((w5(x)-1)-3*w7(x))/2
	carbonexcess=w6(x)-(w7(x)+(w5(x)-1))/2 
	composite=carbonexcess+(0.5*methyleneexcess/25)
	zubarev=w6 - (w7 + (w5-1)) /2
	dowindow/F vankrevH_N
	if (V_FLag!=1)
		execute"vankrevH_N() "
		addmyControls()
		Button button4,pos={390,00},size={90,18},proc=standardHC,title="standard H/C"
	endif
	dowindow/F GrN_C
	if (V_FLag!=1)
		execute"GrN_C()"
		addmyControls()
	endif
	dowindow/F GrH_C
	if (V_FLag!=1)
		execute"GrH_C()"
		addmyControls()
		Button button4,pos={390,00},size={90,18},proc=standardHC,title="standard H/C"	
	endif
	dowindow/F GrDBE
	if (V_FLag!=1)
		execute"GrDBE()"
		addmyControls()
	endif
	dowindow/F tholinomics
	if (V_FLag!=1)
		execute"tholinomics()"
		addmyControls()
	endif
	dowindow/F Grcomposite
	if (V_FLag!=1)
		execute"Grcomposite()"
		addmyControls()
	endif
	dowindow/F GrZubarev
	if (V_FLag!=1)
		execute"GrZubarev()"
		addmyControls()
	endif
	dowindow/F VanKrevO
	if (wavemax(element_o)>0)
		if (V_FLag!=1)
			execute"VanKrevO()"
			addmyControls()
			Button button4,pos={390,00},size={90,18},proc=standardHC,title="standard H/C"	
		endif
	endif
end
function addmycontrols ()
	ControlBar 19
	ModifyGraph cbRGB=(8000,50000,50000)
	Button button1,pos={10,0},size={90,18},proc=removecursorpoint,title="kill this point"
	Button button2,pos={110,0},size={80,18},proc=dorepresentations,title="represent"
	Button button3,pos={200,0},size={80,18},proc=gobackToComp,title="back to comp"
	Button button13,pos={320,00},size={60,18},proc=stopThegeneralHook,title="Hook OFF"

	getwindow /Z kwtopwin title
	if (cmpstr (S_Value, "compare_attributions")!=0)
		SetWindow kwTopWin,hook=get_generalHook,hookevents=1
	endif
end

function standardHC(ctrlName) : ButtonControl
	String ctrlName
	wave hsurc=hsurc
	wave hsurclong=hsurclong
	wave w5=element_H
	wave w6=element_C

	hsurc=(w5(x)-1)/w6(x)
	hsurclong=(w5-1)/w6
	label left "H/C"
end

function gobackToComp(ctrlName) : ButtonControl
	String ctrlName
	dowindow/K vankrevH_N
	dowindow/K vankrevO
	dowindow/K GrN_C
	dowindow/K GrH_C
	dowindow/K GrDBE
	dowindow/K tholinomics
	dowindow/K Grcomposite
	dowindow/K GrZubarev
	startTheHook("bob")
end

function stopThegeneralHook(ctrlName) : ButtonControl
	String ctrlName
	SetWindow kwTopWin, hook=$""
	Button button13,proc=startThegeneralHook,title="Hook ON"
end

function startThegeneralHook(ctrlName) : ButtonControl
	String ctrlName
	SetWindow kwTopWin,hook=get_generalHook,hookevents=1
	Button button13,proc=stopThegeneralHook,title="Hook OFF"
end

Window vankrevH_N() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(1043.25,417.5,1434.75,775.25) hsurc vs nsurc
	ModifyGraph height={Aspect,1},gbRGB=(21760,21760,21760)
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph mrkThick=2
	ModifyGraph zmrkSize(hsurc)={massShort,*,*,8,0.5}
	ModifyGraph zColor(hsurc)={massShort,0,*,YellowHot}
	Label left "(H-N)/C"
	Label bottom "N/C"
	SetAxis/E=1 left 0,2.6
	SetAxis bottom 0,2.6
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=51.90/Y=11.14 trace=hsurc
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1,2,1
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.22820211515864,0.923253150057274
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.48025851938895,0.956013745704467
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.441480611045828,1.67079037800687
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.68954171562867,0.920274914089346
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.24500587544066,1.45635738831615
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.18848413631022,1.49805269186712
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.742420681551117,1.54868270332188
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.543830787309048,1.82565864833906
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.661457109283196,1.71546391752577
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.378848413631022,2.08178694158076

EndMacro


Window GrN_C() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,296,648.75,521) nsurc vs massShort
	ModifyGraph mode=3
	ModifyGraph marker=13
	ModifyGraph msize=4
	ModifyGraph mrkThick=1.5
	Label left "N/C"
	ModifyGraph zColor(nsurc)={DBE,*,*,CyanMagenta}
	SetAxis/A/E=1 bottom
EndMacro

Window GrH_C() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(3,540.5,648,777.5) hsurc vs massShort
	ModifyGraph mode=3
	ModifyGraph marker=44
	ModifyGraph msize=4
	ModifyGraph mrkThick=1.5
	Label left "(H-N)/C"
	ModifyGraph zColor(hsurc)={DBE,*,*,CyanMagenta}
	SetAxis/A/E=1 left
	SetAxis/A/E=1 bottom
EndMacro

Window GrDBE() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,39.5,649.5,276.5) DBE vs massShort
	ModifyGraph mode=3
	ModifyGraph marker=41
	ModifyGraph msize=4
	ModifyGraph mrkThick=1.5
	ModifyGraph zColor(DBE)={DBE,*,*,CyanMagenta}
	SetAxis/A/E=1 bottom
	Label left "DBE=C-H/2+N/2+1"
	TextBox/C/N=text0/F=0/A=MC/X=-44.21/Y=45.33 ""

EndMacro

Window tholinomics() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(1044,39.5,1434,395.75) methyleneexcess vs carbonexcess
	ModifyGraph height={Aspect,1},gbRGB=(21760,21760,21760)
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph mrkThick=2
	ModifyGraph zmrkSize(methyleneexcess)={massShort,*,*,15,1}
	ModifyGraph zColor(methyleneexcess)={massShort,0,*,YellowHot}
	ModifyGraph zero=1
	Label left "methylene excess = (H-3N)/2"
	Label bottom "carbon excess = C-(N+H)/2"
	SetAxis/A/N=1 bottom
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=-41.73/Y=13.49 trace=methyleneexcess
EndMacro

Window VanKrevO() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(1044,39.5,1434,395.75)  HsurClong vs Osurc
	ModifyGraph height={Aspect,1},gbRGB=(21760,21760,21760)
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph mrkThick=2
	ModifyGraph zmrkSize(hsurclong)={massinconnue,*,*,8,0.5}
	Label left "(H-N)/C"
	Label bottom "O/C"
	SetAxis/A/E=1 left
	SetAxis/A/E=1 bottom
	if (wavemax(element_o)>0)
		ModifyGraph zColor(hsurclong)={element_O,*,*,BlueHot}
		ColorScale/C/N=text0/F=0/B=1/A=MC/X=-43.64/Y=-6.78 trace=hsurclong
	endif
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,2,0.204129604559247,0.56232429042305
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.251804610887229,0.484848499298096
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.313916414906079,0.422368022584424
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.358905505384598,0.554826633217408
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.427732099027107,0.714776653604409
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.546919614847062,0.907216521882519
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.580493562965359,1.41705721186608
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.582507999852457,2.00187447390606
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.5,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.559677715132015,1.67197755685787
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.572435815416968,1.15213999060011
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.438475762424962,2.18681668497852
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.371999345150734,2.22930340914382
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.276985071975952,2.20930965659545
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.217894923287749,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.182306538282354,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.156790337712448,2.18181824684143
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.138996145209751,2.21930653286963
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.121201952707053,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.185328193613001,0.512339909052111
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.158469035118363,0.412371146310236
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.137653187285019,0.339893793322376
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.122880650112968,0.274914097540157
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.109115331384466,0.247422687786142
EndMacro

Window Grcomposite() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(654,38.75,1038.75,389.75) composite vs massShort
	ModifyGraph height={Aspect,1},gbRGB=(21760,21760,21760)
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph mrkThick=2
	ModifyGraph zColor(composite)={DBE,*,*,CyanMagenta}
	ModifyGraph zero=1
	ModifyGraph lblMargin(left)=15
	Label left "composite "
	Label bottom "m/z"
	SetAxis/A/E=1 bottom
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=-42.53/Y=13.16 trace=composite
EndMacro

Window GrZubarev() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(653.25,417.5,1037.25,777.5) zubarev vs MassInconnue
	ModifyGraph gbRGB=(52224,52224,52224)
	ModifyGraph mode=3
	ModifyGraph marker=43
	ModifyGraph mrkThick=1.5
	ModifyGraph zero(left)=1
	Label left "z = C - (N + H)/2"
	SetAxis/A/E=1 bottom
	if (wavemax(element_o)>0)
		ModifyGraph zColor(zubarev)={element_O,*,*,BlueHot}
		ColorScale/C/N=text0/F=0/B=1/A=MC/X=-43.64/Y=-6.78 trace=zubarev
	endif
EndMacro

Window table_representations() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(42.75,125.75,1401.75,649.25) qualite,delta,numH,numC,numN,numO,masse,inconnue
	AppendToTable hsurc,nsurc,methyleneexcess,carbonexcess,composite,zubarev
	ModifyTable format(Point)=1,width(numH)=77,width(numC)=59,width(numN)=60,width(numO)=53
EndMacro

//////////////////////////////////////////////////////////////////////////////////////////
// Djevahirdjian Leo 2013                                                     //
//                                                                                      //
// Traitement des spectres, detection des pics, nettoyage    //
// des radios et des arctefacts de ringing.                            //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////


// Réécriture de la procédure développée précédement, prenant en compte les ajouts d'Arthur Choplin.
// Ajout de commentaires, nettoyage du code et réécriture de certains points visant plus d'efficacité.
// Modification de la condition d'Union, utilisation d'un ou exclussif.

// Cette procédure vise à detecter les pics par une structure Union-Find, elle trie ensuite ceux détectés pour
// éliminer ceux issus du bruit radio.

// Fonction de recherche de pics se basant sur Union Find, notée U-F pour la suite, placée juste après cette fonction.
// Une fois les pics détecté un filtre anti-radio est appliqué ce traitement peut être suivit enfin d'un filtre anti-ringing
// codé dans une fonction séparée placée après la structure U-F.
Function detectionUF_et_filtreRadio(masse,intensite, parametreUnion, parametreFiltre, s)
	wave masse,intensite
	string s //set data or roi pour copier le resultat dans l'un ou l'autre
	Variable parametreUnion 	// Parametre mis en place par Arthur Choplin, pilotant clairement les unions quand on se trouve sur un minima. Valeur comprise entre 0 et 1
	Variable parametreFiltre 	// Parametre limite de conservation d'un maxima lors du filtrage. Ne seront conservé que les pics dont le nombre de point nécessaire
							// pour atteindre la moitié de leur intensité est inferieure à ce parametre.
	
	Print "START"
	Print "Demarage de la detection de pics avec la structure Union-Find avec le paramètre : ", parametreUnion
	Print "Detection suivie du filtre anti-radio avec le paramètre : ", parametreFiltre
	Variable debut = DateTime
	
	// Creation de l'index permettant le classement en intensité décroissante,
	// pour parcourir par la suite les données dans un ordre d'intensités décroissantes.
	Make/D/O/N=(numpnts(masse)) indexIntensite
	MakeIndex/R intensite, indexIntensite
	
	// Creation de la Wave contenant les pointeurs (plus précisément leur index) vers les parents du point d'index i.
	// On intitialise cette wave avec des valeurs à -1 pour signifier que le point n'a pas été traité par la strucutre U-F
	Make/D/O/N=(numpnts(masse)) parent
	Variable j
	For(j=0 ; j<numpnts(masse) ; j+=1)
		parent[j] = -1
	endFor
	
	For(j=0 ; j<numpnts(masse) ; j+=1)
	// On initialise le pointeurParent vers le point lui-meme (ce qui est caracteristique d'un maximum, il est le membre de tete du cluster).
	// La suite du code va vérifier que ce point est bien une tete de cluster et si non il réattribuear la valeur du pointeurParent en conséquence.
		parent[indexIntensite[j] ] = indexIntensite[j]
		
		If((intensite[indexIntensite[j] -1] > intensite[indexIntensite[j] ])  && (intensite[indexIntensite[j] +1] > intensite[indexIntensite[j] ]))
			// j est un minima, on va regarder si les clusters à sa droite et à sa gauche doivent fusionner entre eux
			// On commence par chercher les parents du cluster de gauche (clusterA) et du cluster de droite (clusterB).
			// Il s'agit des maximums aux quels sont rattachés les point voisin du minima et tout le cluster associé.
			// Il ne s'agit pas des maxima les plus proches, mais de ceux qui sont des vrais pics et qui ont "gagné" toutes leurs précédentes UNIONs.
			
				Variable parentClusterA = Find(indexIntensite[j] -1,parent)
				Variable parentClusterB = Find(indexIntensite[j]+1,parent)
				
				// CONDITION D'UNION
				// Si pour au moins un des deux pics le rapport Imin/Imax est superieur au parametre d'union ce la signifie qu'il est trop proche du minima,
				// ce n'est pas un pic on doit proceder à une UNION. Avec Imin l'intensité du minima que l'on est en train d'étudier et Imax l'intensité du parent
				// du cluster associé.
				If( ((intensite[indexIntensite[j] ] / intensite[parentClusterA]) > parametreUnion) %^ ((intensite[indexIntensite[j] ] / intensite[parentClusterB]) > parametreUnion) )
					// On rattache le point de minima à un des deux cluster (lequel cela n'a peu d'importance, si l'on veut l'attacher au plus grand cela est possible
					// mais demande un test suplémentaire peu utile qui demande une action de calcul à l'ordinateur. Le test a été fait, on detecte exactement les
					// memes pics que si l'on rattache à celui de plus grande intensité (car ce minima n'a pas d'importance pour la suite de la fonction il n'a plus aucun
					// role dans l'U-F.
					parent[indexIntensite[j] ] = parentClusterA
										
					// On procede à l'union des pics
					union(parentClusterA, parentClusterB)
				Else
					//On a bien deux pics, on rattache le point minima a l'un des deux.
					parent[indexIntensite[j] ] = parentClusterA
				endIf
				
		ElseIf( (intensite[indexIntensite[j] -1] < intensite[indexIntensite[j] ])  && (intensite[indexIntensite[j] +1] < intensite[indexIntensite[j] ]) )
			// Le point est un maxima, on ne fait donc rien et son pointeur parent reste bouclé sur le point lui même
		Else
			// Le point n'est ni un maxima ni un minima ; point intermédiaire dans une pente.
			If(intensite[indexIntensite[j] -1] > intensite[indexIntensite[j] ])
				//Portion decroissante on rattache notre point à gauche
				parent[indexIntensite[j] ] = indexIntensite[j]-1
			Elseif(intensite[IndexIntensite[j]-1] < intensite[indexIntensite[j]])
				//Portion croissante on rattache notre point à droite
				parent[indexIntensite[j] ] = indexIntensite[j]+1
			Else
				// Valeurs avec meme intensite
				// A REFLECHIR
			endIf
		endIf
		
		//Le points extremes du spectres sont conservé et intialisé avec un parent pointant sur eux mêmes.
		parent[0] = 0
		parent[numpnts(masse)-1] = numpnts(masse)-1
		
	endFor
	
	// La fonction d'U-F est terminée elle a detecté les pics, ceux-ci ont un pointeur dans la wave qui parent qui pointe vers eux-même
	// On peut lancer maintenant la récupération de ceux-ci dans une wave de résultat, on couple cette récupération au filtrage anti-bruit radio
	// pour gagner en efficacite.
	
	// Début de la récupération des pics et du filtrage Radio
	
	// Creation des differentes waves necessaire pour la recuperation des pics
	Make/D/O/N=0 masseResult
	Make/D/O/N=0 intensiteResult
	
	// Compteur du nombre de pics detectés par Union-Find avant filtrage.
	// Purement informatif pour l'affichage de la ligne d'information post traitement.
	Variable counterPicsUF
	
	For(j=0 ; j<numpnts(intensite) ; j+=1)
	// On va analyser chaque point considéré comme pic au travers du filtre anti radio. Est considéré comme pic tout point dont le pointeur parent
	// pointe sur le point lui-meme.
	
		If(parent[j] == j)
		// Maximum detecter par la structure Union-Find
		// On demare le filtre qui va compter le nombre de proche voisins avec une intensite plus grande que la moitie de l'intensite du pic.
		// Sera considere comme bruit radio tout pic dont ce nombre est plus grand que la limite fixée par l'utilisateur, car correspondant a un
		// signal avec une forme etalee et non pointue.
			Variable counterVoisins = 0
			Variable i
			// On commence le décompte vers la droite.
			// Le decompte s'arrete quand on passe en dessous de la limite d'intensite.
			// Le decompte s'arette quand le compteur de voisin atteind la limite du filtre car une fois cette limite le pic sera considéré comme du bruit radio.
			// Le decompte s'arette quand on atteind le bout du spectre.
			For(i=0 ; (intensite[j + i] > (intensite[j] / 2)) && (counterVoisins < parametreFiltre) && ((j+i) < numpnts(masse)) ; i+=1)
				counterVoisins += 1
			endFor
			// On commence le decompte vers la gauche, meme logique que precedement
			For(i=0 ; (intensite[j + i] > (intensite[j] / 2)) && (counterVoisins < parametreFiltre) && ((j+i) > 0) ; i-=1)
				counterVoisins += 1
			endFor
			
			If(counterVoisins < parametreFiltre)
			// Si apres les decomptes le pics a toujours un nombre de voisin inferieur a la limite alors c'est un vrai pic et on l'enregistre.
				insertPoints 0,1,masseResult, intensiteResult
				masseResult[0] = masse[j]
				intensiteResult[0] = intensite[j]
			endIf
		
		// Incremente le nombre de pics detectés par l'U-F de 1.
		// Compteur purement informatif pour l'affichage de la ligne d'information post traitement.
		counterPicsUF += 1
		endIf
	endFor
	
	Variable fin = DateTime
	
	// Affichage des informations post traitement
	Print "Recupération des pics et filtrage terminés"
	Print "Nombre de pics détectés par l'U-F avant filtrage anti-radio:", counterPicsUF
	Print "Nombre de pics retenus après filtrage anti-radio :", numpnts(masseResult)
	Print  "\tTemps d'execution : ", (fin-debut)
	Print "END"
	
	// Sauvegarde des waves dans celles d'Atributor et nettoyage de celles utilisées
	Make/D/O masseOld
	Make/D/O intensiteOld
	
	Duplicate/O wave0stripped masseOld
	Duplicate/O wave1stripped intensiteOld
	
	strswitch(s)
		case "ROI":
				Duplicate/O masseResult roi0
				Duplicate/O intensiteResult roi1
				Duplicate/O masseResult roi2
				sort roi0, roi0, roi1
				roi2=defmass(roi0)
			break
		case "data":
				Duplicate/O masseResult wave0stripped
				Duplicate/O intensiteResult wave1stripped	
			break
		default:
			print "No data nor ROI to copy results in."
	endswitch
	
	
	
	KillWaves masseResult intensiteResult

End

	
//////////////////////////////////////////////////////////////////////////
//                    Structure Union-Find                    //
//////////////////////////////////////////////////////////////////////////
	
// Fonction Find, permet de trouver le parent du cluster associé à un point a.
// Lors de la recherche on en profite pour réalocaliser les parents tout au long du parcours
// pour aplatir le chemin d'arborescence du cluster.
Function Find(k,parent)
Variable k
Wave parent
		
	If(parent[k]!=k)
		parent[k] = Find(parent[k],parent)
	endif
		
	return parent[k]
end
	
// Fonction Union, elle fusionne les cluster A et B ensemble en un seul cluster.
// Pour ce faire elle fait pointer le parent du cluster le moins intense vers le parent du cluster
// de plus grande intensité.
Function union(a, b)
	Variable a  // pointeur vers le parent du cluster A
	Variable b  // pointeur vers le parnet du cluster B
	Wave intensite=wave1stripped
	Wave parent=parent
		
	If(intensite[a] > intensite[b])
		parent[b] = a
	Else
		parent[a] = b
	endIf
End

//////////////////////////////////////////////////////////////////////////
//                END Structure Union-Find                 //
//////////////////////////////////////////////////////////////////////////

// Fonction de nettoyage des artefacts de ringing basées sur un algorhytme débile de recherche de symétrie.
// Cette fonction est a appeler après une detection de pic avec Union-Find pour prendre tout son sens.
Function nettoyageRinging(toleranceppm, tailleZoneInspection, valeurConditionIntensiteHaute, valeurConditionIntensiteBasse, limiteConditionIntensite)
	Variable toleranceppm
	Variable tailleZoneInspection
	Variable valeurConditionIntensiteHaute
	Variable valeurConditionIntensiteBasse
	Variable limiteConditionIntensite

	Print "START"
	Print "Démarage du filtre anti-Ringing avec les paramètres suivants :"
	Print "\tTolérance ppm = ", toleranceppm, "\tTaille de la zone d'inspection (uma) = ", tailleZoneInspection
	Print "\tCondition intensité :\tvaleur haute intensité = ", valeurConditionIntensiteHaute, " valeur basse intensité = ", valeurConditionIntensiteBasse
	Print "\t\t\t\t\t\tlimite zones intensité = ", limiteConditionIntensite
	Variable debut = DateTime

	Wave masse = wave0stripped
	Wave intensite = wave1stripped	
	
	Variable nbPicInit = numpnts(masse)
	
	If(masse[0] > masse[1])
		Reverse masse
		Reverse intensite
	endIf
	
	// Creation de l'index permettant le classement en intensité décroissante,
	// pour parcourir par la suite les données dans un ordre d'intensités décroissantes.
	Make/O/N=(numpnts(masse)) indexIntensite
	MakeIndex/R intensite, indexIntensite
	
	// Commence le nettoyage bete et méchant, basé sur une recherche en symétrie classique.
	// On crée la wave qui va nous servir à identifier les Ringings, la seconde va contenir les pointeurs vers les paires.
	Make/O/N=(numpnts(masse)) idRinging, paireRinging
	
	Variable l
	For(l=0 ; l<numpnts(idRinging) ; l+=1)
		idRinging[l] = -1
		paireRinging[l] = -1
	endFor
	
	Variable i
	Variable j
	For(i=0 ; intensite[indexIntensite[i] ] > 0.005 ; i+=1)
		//On inspecte jusqu'a une masse autour du pic central
		For(j=1 ; masse[indexIntensite[i] +j] < (masse[indexIntensite[i] ]+ tailleZoneInspection) ; j+=1)
			Variable masseRecherche = 2*masse[indexIntensite[i] ] - masse[indexIntensite[i] + j ]
			// On recherche la masse attendue, symétrique du point avec une tolérance sur la valeur trouvée. On utilise une function perso dont les résultats sont
			// plus satisfaisants que ceux obtenus avec la recherche d'Igor.
			Variable indexTrouve
			indexTrouve =  trouverValeur(masse, masseRecherche, indexIntensite[i], toleranceppm*masseRecherche/1000000, 1)
			
			If(indexTrouve !=-1)
				// On enregistre les caracteristiques de ces pics symétriques pour les analyser ensuite	
				// Defini une variable memoire
				Variable rapportIntensitePaire = min(intensite[indexTrouve], intensite[indexIntensite[i]+j ])/max(intensite[indexTrouve], intensite[indexIntensite[i]+j ])
				
					If( (intensite[indexTrouve] < intensite[indexIntensite[i]]) && (intensite[indexTrouve] < intensite[indexIntensite[i]]) )
					// Ce test evite les morts à la con, les repliques doivent etre plus petite que le pic central
					
						If( (idRinging[indexTrouve] == -1) && (idRinging[indexIntensite[i]+j ] ==-1) ) // les points n'ont pas déja de paire
						
							idRinging[indexTrouve] =rapportIntensitePaire
							paireRinging[indexTrouve] = indexIntensite[i]+j
							
							idRinging[indexIntensite[i]+j ] = rapportIntensitePaire
							paireRinging[indexIntensite[i] + j ] = indexTrouve
							
						Else // un des points appartient déjà à une paire il faut vérifier si on a pas trouvé une meilleure paire
						
							If(idRinging[indexTrouve] != -1)
								If(idRinging[indexIntensite[i]+j ] != -1) // Les deux ont déjà une paire
								
									If( (rapportIntensitePaire > idRinging[indexTrouve]) && (rapportIntensitePaire > idRinging[indexIntensite[i]+j ]) )
										// Plus interessant pour les deux pics, ils sont mieux ensemble qu'avec leur ancien colègue (rapport intensité plus grand)
										idRinging[indexTrouve] = rapportIntensitePaire
										paireRinging[indexTrouve] = indexIntensite[i]+j
										
										idRinging[indexIntensite[i]+j ] = rapportIntensitePaire
										paireRinging[indexIntensite[i] + j ] = indexTrouve
										
										// Annule la précédante detection sur les pics associés
										idRinging[paireRinging[indexTrouve]] = -1
										paireRinging[paireRinging[indexTrouve]] = -1
										idRinging[paireRinging[indexIntensite[i]+j]] = -1
										paireRinging[paireRinging[indexIntensite[i]+j]] = -1
										
										Else
										//Plus interessant pour seulement un des deux, pour le moment on ne fait rien
									endIf
								
								Else //c'est le pic en indexTrouve qui a deja une paire, et lui seul
									If(rapportIntensitePaire > idRinging[indexTrouve])
									
										// Annule la précédante detection sur le pic associé
										idRinging[paireRinging[indexTrouve]] = -1
										paireRinging[paireRinging[indexTrouve]] = -1
									
										idRinging[indexTrouve] = rapportIntensitePaire
										paireRinging[indexTrouve] = indexIntensite[i]+j
										
										idRinging[indexIntensite[i]+j ] = rapportIntensitePaire
										paireRinging[indexIntensite[i] + j ] = indexTrouve										
										
									endIf
								endIf
								
							Else // c'est le pic en indexInt[i] + j qui a deja une paire, et lui seul
								If(rapportIntensitePaire > idRinging[indexIntensite[i] + j ])
								
								 	// Annule la précédante detection sur le pic associé
									idRinging[paireRinging[indexIntensite[i] + j ]] = -1
									paireRinging[paireRinging[indexIntensite[i] + j ]] = -1
																						
									idRinging[indexTrouve] = rapportIntensitePaire
									paireRinging[indexTrouve] = indexIntensite[i]+j
										
									idRinging[indexIntensite[i]+j ] = rapportIntensitePaire
									paireRinging[indexIntensite[i] + j ] = indexTrouve
									
								endIf
							endIf
						
						endIf
					endIf
			Else
				// pas de symétrique ce point n'est pas ringing		
			endIf
		endFor
	endFor
	
	// On commence à filtrer pour recuperer les pics non ringing	
	Make/O/N=0 massePicNonRinging
	Make/O/N=0 intensitePicNonRinging
	
	Variable counter=0
	For(i=0 ; i<numpnts(idRinging) ; i+=1)
		If(idRinging[i] == -1)
			insertPoints 0, 1, massePicNonRinging, intensitePicNonRinging
			
			massePicNonRinging[0] = masse[i]
			intensitePicNonRinging[0] = intensite[i]
		Else
			If(intensite[i]>limiteConditionIntensite)
				// Pic symétrique intense
				If(idRinging[i]<valeurConditionIntensiteHaute)
					// Le rapport en intensité du couple symétrique est trop faible (pic pas assez corrélés en intensité),
					// ce n'est donc pas du ringing et on le récupère.
					insertPoints 0, 1, massePicNonRinging, intensitePicNonRinging
			
					massePicNonRinging[0] = masse[i]
					intensitePicNonRinging[0] = intensite[i]
				endIf		
			Else
				// Pic symétrique peu intense
				If(idRinging[i]<valeurConditionIntensiteBasse)
					// Le rapport en intensité du couple symétrique est trop faible (pic pas assez corrélés en intensité),
					// ce n'est donc pas du ringing et on le récupère. On utilise la conditon en intensité souple car le pic
					// est peu intense donc sujet au bruit.
					insertPoints 0, 1, massePicNonRinging, intensitePicNonRinging
			
					massePicNonRinging[0] = masse[i]
					intensitePicNonRinging[0] = intensite[i]
				endIf
			endIf
		endIf
	endFor
		
	Variable fin = DateTime
	
	Print "Filtrage terminé."
	Print "\tNombre de pics ringing detectés : ", (nbPicInit - numpnts(massePicNonRinging)), "\tTemps d'execution : ", (fin-debut)
	Print "END"
	
	// Copie les résultats dans les waves d'attributor
	Duplicate/O massePicNonRinging wave0stripped
	Duplicate/O intensitePicNonRinging wave1stripped
		
	// Fait le menage
	Killwaves massePicNonRinging, intensitePicNonRinging, idRinging, paireRinging
End


// Cherche le point le plus proche de la valeur recherchée, cette fonction vient remplacer la fonction Igor qui n'était pas satisfaisante.
// Retroune l'index dela valeur trouvée si celle-ci est dans l'intervalle de tolérance, sinon renvoit la valeur -1.
Function trouverValeur(waveSource, valeurRecherche, indexDepart, tolerance, boolReverseMode)

	Wave waveSource
	Variable valeurRecherche
	Variable indexDepart // Position de départ dans la wave source pour la recherche
	Variable tolerance // Tolerance sur la masse trouvée
	Variable boolReverseMode // Prend la valeur 1 si l'on travaille dans le sens décroissant.
	
	Variable ecart1 = abs(valeurRecherche - waveSource[indexDepart]) // Contient l'ecart de la valeur en cours d'étude
	Variable ecart2 = 1000000
	
	// On défini le signe de l'incrément pour la boucle de recherche (ce signe pilote le sens de recherche)
	Variable increment
	If(boolReverseMode == 1)
		increment = -1
	Else
		increment = 1
	endIf
	
	Variable mem = -5
	
	// On commence la recherche, tant que l'écart du nouveau point est plus petit que celui precedent on continue la recherche.
	// Des que l'on trouve que ce n'est plus le cas (on a depassé la valeur la plus proche) on arrete de chercher pour passer à la suite
	Variable i
	For(i= increment ; ecart2>ecart1 ; i+=increment)
		ecart2 = ecart1
		ecart1 = abs(valeurRecherche - waveSource[indexDepart+i])	
		// On enregistre la position du point précédent (qui sera la valeur de retour que l'on veut en sortie, car il ne faut pas oublier que
		// l'on sort de la boucle un tour plus tard (on fait l'enregistrement dans la boucle for car en dehors la variable i peut être détruite).	
		mem = indexDepart + i - increment
	endFor
	
	// On a trouvé le plus proche on regarde si il est dans l'intervale de tolérance
	Variable test = abs(waveSource[mem] - valeurRecherche)
	If(test < tolerance)
		// Valeur dans l'intervalle, elle est donc acceptable
		return mem
	Else
		// Valeur non acceptable on retourne -1
		return -1
	endIf
End


//////////////////////////////////////////////////////////////////////////
//                        GUI personelle                        //
//////////////////////////////////////////////////////////////////////////
								
										
//Function creant le panel ou apparaisse les boutons dont j'ai besoin				

Function filtresGUI()
	dowindow/K Filtres
	NewPanel/W=(952,707,1152,784) as "Filtres"
	Button Bouton1, proc=bouton1, pos={5,5}, size={190,30}, title="Union Find + radio"
	Button Bouton2, proc=bouton2, pos={5,40}, size={190,30}, title="anti-Ringing"
End

// Liste des boutons de controle

Function bouton2(NomButton) : ButtonControl
	String NomButton
	
	nettoyageRinging(2,0.3,0.8,0.2,0.00005)
End

//////////////////////////////////////////////////////////////////////////
//                    END GUI personelle                     //
//////////////////////////////////////////////////////////////////////////


Function ButtonPeaks(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave/Z wave2stripped, wave0stripped
	switch( ba.eventCode )
		case 2: // mouse up
			//keepeaks()   						//mis en commentaire par Roland en Aout 2013 pour lancer plutot le nouveau traitement ALAL
			//wave2stripped=defmass(wave0stripped)
			filtresGUI()
			break
	endswitch

	return 0
End

/////////////////////////////////////////////////////////////////////////////////
///////Corrections de Roland le 6 novembre 2013/////////////
/////////////////////////////////////////////////////////////////////////////////
function gicle2data(newname)
string newname
string les1, les0
//wave intensitegiclee intensitegiclee
//wave massegiclee massegiclee
les1=newname+"_1stripped"
les0=newname+"_0stripped"
duplicate/O intensitegiclee $les1
duplicate/O massegiclee $les0
sort $les0,  $les0, $les1
listelesdonnees()
end

function goods2data(newname)
string newname
string les1, les0
les1=newname+"_1stripped"
les0=newname+"_0stripped"
duplicate/O Intensitepic $les1
duplicate/O massinconnue $les0
sort $les0,  $les0, $les1
listelesdonnees()
end
///////////////////////////////////////////////////////////////////////////////////


///////////////correction introduisant les centredemasses pour recalibration 7 novembre 2013

threadsafe function generelesmat(kprem)
variable kprem
variable t=ticks
wave/Z wave0stripped, wave2stripped, dataindex, wave1stripped
Make/FREE/O/N=(kprem) dim=0
MatrixOP/FREE/O cropmassx=dim x dim^t
duplicate/O cropmassx, diffm, diffdm, diffint, centremasse, slopes, modules, intercepts//, calcerror, adja//génère la matrice de tous les écart de defaut de masse point à point, paramètre de sortie
Duplicate/FREE/O cropmassx, cropdefmassx, cropintx
Multithread cropmassx=wave0stripped[dataindex[x]]
Multithread cropdefmassx=wave2stripped[dataindex[x]]
Multithread cropintx=wave1stripped[dataindex[x]]
MatrixOP/FREE/O cropmassy=cropmassx^t
MatrixOP/FREE/O cropdefmassy=cropdefmassx^t
MatrixOP/FREE/O cropinty=cropintx^t
FastOP diffm=cropmassx-cropmassy //calcul
FastOP diffdm=cropdefmassx-cropdefmassy//calcul
FastOP diffint=cropintx+cropinty
FastOP centremasse=0.5*cropmassx+0.5*cropmassy
Fastop slopes = diffdm/diffm //calcul des pentes
FastOP intercepts=cropdefmassy-slopes*cropmassy
FastOP modules=diffm*diffm
FastOP modules=modules+diffdm*diffdm
Multithread modules = sqrt(modules)
//make/O/B/N=(kprem,kprem,114) slopetrix
//Multithread slopetrix=0
//make/O/T/N=(kprem,kprem) slopenames
//slopenames=""
print "generelesmat:", ticks-t
end

function getslopestat(slo)
variable slo
variable k,n,t=ticks
wave/Z slopes, modules, wave0stripped,dataindex,wave1stripped,intercepts, centremasse
make/O histoslo // génère les histogramme associés, paramètres de sorti
histogram/B={-0.01,2e-5,1000} slopes histoslo
duplicate/O slopes slopematch
slopematch=(slopes>slo && slopes<=slo+2e-5)
Matrixop/O checkslo=(sumcols(slopematch))^t
//simplifysmatch()
Make/O/N=(numpnts(modules)) modtrans, pentrans, interceptrans, centremassetrans, xxtrans
pentrans=slopematch[x]*slopes[x]
modtrans=slopematch[x]*modules[x]
interceptrans=slopematch[x]*intercepts[x]
centremassetrans=slopematch[x]*centremasse[x]
xxtrans=slopematch[x]*x
sort modtrans, modtrans, pentrans, interceptrans, xxtrans, centremassetrans
deletepoints 0,numpnts(modules)-sum(checkslo), modtrans, pentrans, interceptrans, xxtrans, centremassetrans
//sort/R modtrans, modtrans, pentrans, interceptrans
Make/O/N=0 roi0,roi1,roi2
n=numpnts(checkslo)
	for(k=0;k<n;k+=1)
		if(checkslo[k]!=0)
		insertpoints 0,1,roi0,roi1,roi2
		roi0[0]=wave0stripped[dataindex[k]]
		roi1[0]=wave1stripped[dataindex[k]]
		endif
	endfor
	sort roi0, roi0,roi1
	roi2=roi0-round(roi0)
end
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////Calibration non linéaire à 2 points
////////////////////////////////////////////////////////////////////////
function delM(M, calwidth)//interpoler delM, l'integrer et corriger avec 2 points.
variable M, calwidth
wave diffm, centremasse, wave0stripped, w_coef, pointconf
variable reso, B,A, errmes1, errmes2, errcalc1, errcalc2
nvar degreedelm
duplicate/O wave0stripped error0 error1
if(M>0)
make/O histo
histogram/B={M-0.5,0.5,2} diffm histo
reso=sum(histo)
histogram/B={M-0.5,1/reso,reso} diffm histo
findlevel histo wavemax(histo)
M=V_levelX
histogram/B={M-calwidth*wavemax(histo)/reso,2*calwidth/reso,wavemax(histo)} diffm histo
duplicate/O diffm diffmmatch
diffmmatch = (diffm>(M-calwidth*wavemax(histo)/reso) && diffm<(M+calwidth*wavemax(histo)/reso))
make/O/n=(numpnts(diffm)), diffmcrop, centremasscrop
diffmcrop=diffmmatch[x]*diffm[x]
centremasscrop=diffmmatch[x]*centremasse[x]
sort centremasscrop, centremasscrop, diffmcrop
deletepoints 0,numpnts(diffm)-sum(diffmmatch), diffmcrop, centremasscrop
duplicate/O diffmcrop diffmcropmask
CurveFit/Q/M=2/W=0 poly (degreedelm), diffmcrop/X=centremasscrop/D/M=diffmcropmask
//error1=(W_coef[3]*error0^4/4+W_coef[2]*error0^3/3+W_coef[1]*error0^2/2)/V_levelx
error1=poly(W_coef,error0[x])
differentiate error1 /X=error0 /D=error1
integrate/METH=1 error1 /X=error0 /D=error1
integrate/METH=1 error1 /X=error0 /D=error1
error1/=V_levelx
else
error1=0
endif
errmes1=0+pointconf[0][0]-pointconf[1][0]
errmes2=0+pointconf[0][1]-pointconf[1][1]
findvalue/V=(pointconf[0][0]) error0
errcalc1=error1[V_value]
findvalue/V=(pointconf[0][1]) error0
errcalc2=error1[V_value]
B=(errmes2-errcalc2-errmes1+errcalc1)/(pointconf[0][1]-pointconf[0][0])
A=errmes1-errcalc1-B*pointconf[0][0]
print B,A
error1=error1+B*error0+A
killwaves diffmmatch
end

///révision le 28 janvier de valse pour inclure le calcul de la lorentzienne
function valse(kprem) // reviewed 13 avril 2012, segmentation en generelesmat()
variable kprem
variable t=ticks
wave/Z wave0stripped, wave2stripped, dataindex, wave1stripped
	if (kprem<1)
		kprem=numpnts(wave0stripped)
	endif
controlinfo/W=advancedmanager TabCon_List_DATA //récupère le nom de l'échantillon
string laun, ladeux, sxplode, syplode
wave/T listwave=$S_value
variable taille, contrainte, maxx,maxy, n=1000
laun=listwave[V_value][1]
ladeux="Histoiso_"+laun
sxplode="Xplode_"+laun
syplode="Yplode_"+laun
laun="Histoslo_"+laun
generelesmat(kprem)
make/O histoslo, histomo, histoiso // génère les histogramme associés, paramètres de sorti
//histogram/B={-0.01,0.02/n,n} slopes histoslo //calcul l'histogramme des pentes
//histogram/B={1,1,100} modules histomo //calcul l'histogramme des modules
histogram/B={.5,2/100000,100000} diffm histoiso
print ticks-t,"valse-histogrammes";t=ticks
//killwaves cropmass, cropdefmass, diffm,diffdm,slopes,modules,diffint // détruit des grosses wave/Z de calcul
variable v_fiterror=0
CurveFit/Q/M=2/W=2 lor, histoslo/D
wave W_coef
print ticks-t,"valse-fitCauchy";t=ticks
make/O/N=(n) cauchyx, cauchyy, cauchyxplode, cauchyyplode
cauchyx=pnt2x(histoslo,p)
cauchyy=(W_coef[0]+W_coef[1]/((cauchyx-W_coef[2])^2+W_coef[3]))
duplicate/O histoslo xplode yplode, rawxplode, rawyplode
rawxplode=(histoslo)*cos(atan(pnt2x(histoslo,p)))
rawyplode=(histoslo)*sin(atan(pnt2x(histoslo,p)))
cauchyxplode=cauchyy*cos(atan(pnt2x(histoslo,p)))
cauchyyplode=cauchyy*sin(atan(pnt2x(histoslo,p)))
xplode=max(0,(histoslo-cauchyy))*cos(atan(pnt2x(histoslo,p)))
yplode=max(0,(histoslo-cauchyy))*sin(atan(pnt2x(histoslo,p)))
maxx=max(wavemax(xplode), abs(wavemin(xplode)))
maxy=max(wavemax(yplode), abs(wavemin(yplode)))
contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
xplode*=contrainte
yplode*=contrainte
cauchyxplode*=contrainte
cauchyyplode*=contrainte
rawxplode*=contrainte
rawyplode*=contrainte
xplode+=wavemin(wave0stripped)
cauchyxplode+=wavemin(wave0stripped)
rawxplode+=wavemin(wave0stripped)
setaxis/W=dmvm/A bottom
print ticks-t,"valser";t=ticks
end

//revue le 28 janvier pour adapté à la révision de valse qui soustray la Cauchy law.
function getdmvmslopemode(info)
string info
wave/Z yplode,xplode, incoxdm, incox, histoslo
make/o/d/n=1 pozmassedmvm, pozdmassedmvm
variable slope, dabin
pozmassedmvm=axisvalfrompixel("dmvm","bottom",str2num(stringbykey("mousex",info)))
pozdmassedmvm=axisvalfrompixel("dmvm","left",str2num(stringbykey("mousey",info)))
duplicate/O xplode distplode
distplode=(((pozdmassedmvm-yplode))^2+(0.1*(pozmassedmvm-xplode))^2)
findlevel/Q/P distplode wavemin(distplode)
incoxdm=yplode[v_levelx]
incox=xplode[v_levelx]
slope=pnt2x(histoslo,v_levelx)
doupdate
	if(stringmatch(stringbykey("event",info),"mousedown")==1)
		print slope
		getslopestat(slope)
		findmostzero()
		analysemult(1)
	endif
end


//révision 28 janvier 2014 avec tri en masse obligatoire (réordonne les données par masse croissante)
function activeladata(name)// ATTENTION 100 PARAMETRE ARBITRAIRE
string name
variable t=ticks
nvar inco, intenscumthres
if(cmpstr(name,"")==0)
	print "nom invalide"
	return -1
endif
						
string laun, lazero,latrois
wave/Z topbrighty, topbrightx, incoy, incox, logwave1, incoxdm
laun=name+"_1stripped"
lazero=name+"_0stripped"
latrois=name+"_3stripped"
Duplicate/O $laun wave1stripped, dataindex
Duplicate/O $lazero wave0stripped
Duplicate/O $lazero wave2stripped
if(waveexists($latrois))
	Duplicate/O $latrois wave3stripped
endif
duplicate/O wave1stripped logwave1
Sort wave0stripped, wave0stripped, wave1stripped,wave3stripped
multithread wave2stripped=defmass(wave0stripped)
multithread logwave1=log(wave1stripped)
makeindex/R wave1stripped, dataindex
incoy=wave1stripped[dataindex[0]]
incox=wave0stripped[dataindex[0]]
incoxdm=defmass(incox)
inco=incox[0]
Variable twosig
twosig=intenstat(intenscumthres)
//valse(twosig)
generelesmat(twosig)
print ticks-t,"valser";t=ticks
setaxis/W=agregateur bottom wavemin(wave0stripped), wavemax(wave0stripped)
DoUpdate
setaxis/W=agregateur /A=2 left
topbright()
calcinfo()
print ticks-t,"calculer";t=ticks
end 


////introduit le 28 janvier 2014 pour la calibration à 2 points
function mainmass()
wave simu_mass, simu_proba
variable k
duplicate/O simu_proba, transindex
makeindex/R simu_proba, transindex
k=simu_mass[transindex[0]]
killwaves transindex
return k
end


////introduit le 28 janvier 2014 pour la calibration à 2 points
function StoConfPnt()
wave simu_mass
variable/g inco
if (waveexists(pointconf)==0)
make/D/O/N=(2,2) pointconf
pointconf[0][0]=inco ; pointconf[1][0]=mainmass()
else
insertpoints/m=1 0,1, pointconf
pointconf[0][0]=inco ; pointconf[1][0]=mainmass()
deletepoints/M=1 2,1, pointconf
endif
end

Function ButtonProc_2(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			stoconfpnt()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

////introduit le 28 janvier 2014 pour la calibration à 2 points
Window HiResCalMain() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(205.5,60.5,885.75,509) error1 vs error0
	AppendToGraph/L=upleft diffmcrop vs centremasscrop
	AppendToGraph/L=upleft fit_diffmcrop
	AppendToGraph pntserror vs pntsmes
	AppendToGraph fit_pntserror
	ModifyGraph mode(diffmcrop)=3,mode(pntserror)=3
	ModifyGraph marker(pntserror)=59
	ModifyGraph lSize(error1)=2
	ModifyGraph rgb(error1)=(65280,0,52224),rgb(fit_diffmcrop)=(0,0,65280),rgb(pntserror)=(65280,0,0)
	ModifyGraph rgb(fit_pntserror)=(0,0,65280)
	ModifyGraph grid(upleft)=1
	ModifyGraph zero(left)=1
	ModifyGraph axOffset(left)=21.2222
	ModifyGraph gridRGB(upleft)=(0,0,0)
	ModifyGraph lblPos(left)=70,lblPos(upleft)=75
	ModifyGraph freePos(upleft)=0
	ModifyGraph axisEnab(left)={0,0.45}
	ModifyGraph axisEnab(upleft)={0.55,1}
	Label left "Mass Correction (u)"
	Label upleft "Values around mass gap (u)"
	ValDisplay valdisp0,pos={10,86},size={150,14},title="In data",format="%15.13f"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"pointconf[0][0]"
	ValDisplay valdisp1,pos={10,103},size={150,14},title="In theory"
	ValDisplay valdisp1,format="%15.13f",limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"pointconf[1][0]"
	GroupBox group0,pos={5,68},size={165,55},title="Point of confidence 1"
	Button button0,pos={4,7},size={167,49},proc=ButtonProc_2,title="Set current attribution\r as point of confidence 1 \r (passes POC1 to POC 2)"
	Button button0,fColor=(32768,65280,0)
	ValDisplay valdisp2,pos={10,152},size={150,14},title="In data",format="%15.13f"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000},value= #"pointconf[0][1]"
	ValDisplay valdisp3,pos={10,169},size={150,14},title="In theory"
	ValDisplay valdisp3,format="%15.13f",limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp3,value= #"pointconf[1][1]"
	GroupBox group1,pos={5,134},size={165,55},title="Point of confidence 2"
	GroupBox group2,pos={5,193},size={164,58},title="Frequent gap"
	SetVariable setvar0,pos={15,210},size={150,16},proc=SetVarProc_1,title="Mass gap"
	SetVariable setvar0,valueBackColor=(65280,43520,0),value= massgap
	SetVariable setvar1,pos={15,230},size={150,16},proc=SetVarProc_1,title="Gap width (arbitrary)"
	SetVariable setvar1,valueBackColor=(65280,43520,0),value= calwidth
	Button button1,pos={4,326},size={167,49},proc=ButtonProc_4,title="Correct mass"
	Button button1,fColor=(65280,43520,0)
	SetVariable setvar2,pos={5,256},size={218,16},title="Fit Polynomial with degree of liberty:"
	SetVariable setvar2,limits={3,10000,1},value= degreedelm
	//Button button2,pos={4,280},size={167,42},proc=ButtonProc_7,title="Load information from\rthe current molecule list"
	//Button button2,fColor=(32768,65280,0)
	PopupMenu popup1, win=HiResCalMain,pos={4,280},title="Use this list",bodywidth=0,popvalue="Calibrator", value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")",proc=PopMenuProc_HiResCal
	Button button3,pos={4,379},size={167,45},proc=ButtonProc_5,title="ROI 2 data"
	Button button3,fColor=(65280,43520,0)
EndMacro

Function PopMenuProc_HiResCal(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			CalFromMolbag(popStr)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

////introduit le 28 janvier 2014 pour la calibration à 2 points
function initiatecal()
wave wave0stripped, pointconf
variable errmes1, errmes2, errcalc1, errcalc2, A, B
variable/G degreedelm
duplicate/O wave0stripped error0 error1 
make/O/N=0 centremassecrop, diffmcrop, fit_diffmcrop, pntserror, pntsmes
errmes1=0+pointconf[0][0]-pointconf[1][0]
errmes2=0+pointconf[0][1]-pointconf[1][1]
findvalue/V=(pointconf[0][0]) error0
errcalc1=error1[V_value]
findvalue/V=(pointconf[0][1]) error0
errcalc2=error1[V_value]
B=(errmes2-errcalc2-errmes1+errcalc1)/(pointconf[0][1]-pointconf[0][0])
A=errmes1-errcalc1-B*pointconf[0][0]
error1=B*error0+A
end

////introduit le 28 janvier 2014 pour la calibration à 2 points
function correctmass()
wave wave0stripped, wave1stripped, wave2stripped, roi0, roi1, roi2, error1
duplicate/O wave0stripped, roi0
duplicate/O wave1stripped, roi1
duplicate/O wave2stripped, roi2
roi0-=error1
roi2=defmass(roi0)
end

////introduit le 28 janvier 2014 pour la calibration à 2 points
Function ButtonProc_3(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			execute "HiResCalMain()"
			initiatecal()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

////introduit le 28 janvier 2014 pour la calibration à 2 points
Function SetVarProc_1(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			nvar massgap, calwidth
			delM(massgap,calwidth)
		case -1: // control being killed
			break
	endswitch

	return 0
End

////introduit le 28 janvier 2014 pour la calibration à 2 points

////introduit le 28 janvier 2014 pour la transformation des roi par boutons

////introduit le 28 janvier 2014 pour la transformation des roi par boutons
Function ButtonProc_6(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			invertroi()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc_2(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
nvar simu_prof_sigma
simu_prof_sigma=mixsimu()
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Menu "Project"
		"Import Project", loadproject()
		"Save Project", saveproject()
		"Restore default windows layout", setLayaout(windefaultlayout)
		"Layout for small screens", setLayaout(lowreslayout)
		"Save current windows layout", storeLayout()
		"Restore last custom windows layout", setLayaout(WinLayout)
end
 
//revue le 6 fevrier 2014 par Roland pour la bonne gestion des VanK
function removecursorpoint(ctrlName) : ButtonControl
                String ctrlName
                pauseupdate
                variable pt=pcsr(A)
                wave w1=proximite
                wave w2=massinconnue
                wave w3=intensitepic
                wave w4=delta
                wave w5=element_H
                wave w6=element_C
                wave w7=element_N
                wave w8=element_O
                wave gicm=massegiclee
                wave gicint=intensitegiclee
                insertpoints 0,1, gicm, gicint
                gicm[0]=w2[pt]
                gicint[0]=w3[pt]
                deletepoints pt, 1, w1, w2, w3, w4, w5, w6, w7, w8
                getwindow /Z kwtopwin title
                if (cmpstr (S_Value, "compare_attributions")==0)
                               getcursorpoint("bob")
                else
                               dorepresentations ("bob")
                endif
                sort gicm, gicm, gicint
                goods2data("bel");gicle2data("gic");
end

//revue le 6 fevrier 2014 pour zoom scroll//revue le 6 février 2018 pour named hookfunction
function getdmvmMSMSmode(s)
	STRUCT WMWinHookStruct &s
//print s
variable hookres = 0
nvar findindata, treedepth,inco
nvar intervalleHoughSlopes
nvar borneHoughSlopes1
nvar borneHoughSlopes2
nvar intervalleHoughIntercepts
nvar borneHoughIntercepts1
nvar borneHoughIntercepts2
nvar storageHoughResolution
wave/Z wave2stripped, wave1stripped, wave0stripped
wave currentseg, incoxdm, incox, incoy
variable/G pozmassedmvm=s.mouseLoc.h, pozdmassedmvm=s.mouseLoc.v, lastx, lasty, zoomxfactor
pozmassedmvm=axisvalfrompixel("dmvm","bottom",pozmassedmvm)
pozdmassedmvm=axisvalfrompixel("dmvm","left",pozdmassedmvm)
//print info
setwindow dmvm, hookcursor=0
if(s.eventCode==5 && numpnts(roi0)==0)
makeseg()
elseif(s.eventCode==3  || s.eventCode==5)
	if(s.eventMod==16 || s.eventMod==17)
	hookres=1
	popupcontextualmenu "Add current m/z to ROI;Attribute segment;Streak to ROI;Tree to ROI;Calculate"
	strswitch(s_selection)
		case "Add current m/z to ROI":
			AddCurrent2ROI()
			break
		case "Attribute segment":
			attribcombi(sqrt((currentseg[0][1]-currentseg[1][1])^2+(currentseg[0][0]-currentseg[1][0])^2),(currentseg[0][0]-currentseg[1][0])/(currentseg[0][1]-currentseg[1][1]))
			break
		case "Streak to ROI":
			dirstreak2ROI()
			break
		case "Tree to ROI":
			tree2roi()
			break
		case "Calculate":
			nvar charge
			Calculate(inco,incoy[0],charge)
			break
	endswitch
	return hookres
	else
	make/O/D/N=1 modufromdmvm,slofromdmvm, interceptfromdmvm, plotslodmvm, plotinterdmvm
	currentseg[0][0]=incoxdm[0]
	currentseg[0][1]=incox[0]
	currentseg[0][2]=incoy[0]
	findclosepntdmvm(pozmassedmvm,pozdmassedmvm)
	findindata=1
	controlinfo/W=dmvm popup0
	generatetree(S_value,treedepth,inco)
	currentseg[1][0]=incoxdm[0]
	currentseg[1][1]=incox[0]
	currentseg[1][2]=incoy[0]
	modufromdmvm=sqrt((currentseg[1][0]-currentseg[0][0])^2 + (currentseg[1][1]-currentseg[0][1])^2)
	slofromdmvm= ((currentseg[1][0]-currentseg[0][0])/(currentseg[1][1]-currentseg[0][1])) / (abs((currentseg[1][0]-currentseg[0][0])/(currentseg[1][1]-currentseg[0][1]))<max((borneHoughSlopes1),(borneHoughSlopes2)) )
	interceptfromdmvm=currentseg[0][0]-slofromdmvm*currentseg[0][1]
	plotslodmvm=(slofromdmvm-borneHoughSlopes1)*(storageHoughResolution/intervalleHoughSlopes)
	plotinterdmvm=(interceptfromdmvm-borneHoughIntercepts1)*(storageHoughResolution/intervalleHoughIntercepts)
	majMarkerFindOrNot()
	endif
endif
if(s.eventCode==4 && numpnts(cursint)==1)
killcontrol sliderzoom
lastx=s.mouseLoc.h
lasty=s.mouseLoc.v
endif
if(waveexists(currentseg) && currentseg[1][0]!=wave2stripped[updateDistanceDMVM(pozmassedmvm,pozdmassedmvm)] && s.eventCode==4)
make/O/D/N=1 modufromdmvm,slofromdmvm, interceptfromdmvm, plotslodmvm, plotinterdmvm
//attribcombi(modu,slo)
//print modu,slo
elseif(waveexists(currentseg)==0)
make/o/n=(2,3) currentseg
endif
if(s.eventCode==22)
	zoomHook(s.winName,"bottom",pozmassedmvm,0.1*s.wheelDy,wavemin(wave0stripped),wavemax(wave0stripped))
	zoomHook(s.winName,"left",pozdmassedmvm,0.1*s.wheelDy,wavemin(wave2stripped),wavemax(wave2stripped))
	if(s.eventMod==8)
		getaxis/Q/W=$s.winName bottom
		setaxis/W=agregateur bottom V_min, V_max
	endif
	topbright()
	calcinfo()
return 0
endif
return 0
end

////revue le 28 janvier 2014 pour la navigation par zoom à la molette
//revue le 6février2018 pour généralisation et compatibilité avec les hook named
Function zoomHook(win,axe,pos,factor,binf,bsup) //pos en valeur d'axe
string win,axe
variable pos,factor,binf,bsup
variable newmin, newmax
getaxis/Q/W=$win $axe
newmax=V_max
newmin=V_min
newmax=min(pos+(newmax-pos)*(1-factor),bsup)
newmin=max(pos-(pos-newmin)*(1-factor),binf)
setaxis/W=$win $axe newmin, newmax
End

Function SliderZoomProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
variable newleft, newright, newtop, newbottom, pozmasse, pozint
getaxis/Q/W=agregateur bottom
newright=V_max
newleft=V_min
nvar lastx, lasty
wave wave0stripped
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				pozmasse=axisvalfrompixel("agregateur","bottom",lastx)
				newright=min(pozmasse+(newright-pozmasse)*(1+0.2*sign(-sa.curval)),wavemax(wave0stripped))
				newleft=max(pozmasse-(pozmasse-newleft)*(1+0.2*sign(-sa.curval)),wavemin(wave0stripped))
				//print pozmasse, newleft, newright
				setaxis/W=agregateur bottom newleft, newright
				setaxis/W=agregateur/A=2 left
				setaxis/W=dmvm bottom newleft, newright
				setaxis/W=dmvm/A=2 left
				if(itemsinlist(winlist("PeakQualCont",";","")))
					setaxis/W=PeakQualCont bottom newleft, newright
					setaxis/W=PeakQualCont top newleft, newright
				endif
				topbright()
				calcinfo()
			endif
			break
	endswitch

	return 0
End

function addROIfromagreg()
wave roi0, roi1, roi2, wave0stripped, wave1stripped,wave2stripped, roipnts
getmarquee/W=agregateur bottom
variable n=mass2pnt(V_right)-mass2pnt(V_left)
//print v_right,v_left, mass2pnt(V_right),mass2pnt(V_left)
if(stringmatch(S_marqueeWin,"agregateur"))
	insertpoints 0,n,roipnts,roi0, roi1, roi2
	roipnts[0,n-1]=round(mass2pnt(V_left))+x
	roi0[0,n-1]=wave0stripped[roipnts[x]]
	roi1[0,n-1]=wave1stripped[roipnts[x]]
	roi2[0,n-1]=wave2stripped[roipnts[x]]
	sort roi0, roi0, roi1, roi2
	calcinfo()
endif
end

///////////////////Peakfit proto
function pearson4(w,x) : fitfunc
wave w
variable x
return w[0]*(1+((x-w[1])/w[2])^2)^(-w[3])*exp(-w[4]*atan(((x-w[1])/w[2])))
//return w[0]*exp(-0.5*((ln(x)-w[1])/(w[3]))^2)
end

function mixer(ww,xx)
wave ww
variable xx
variable k, n=max(1,dimsize(ww,1)), m=dimsize(ww,0), val=0
make/O/D/N=(m) transcoef
	for(k=0;k<n;k+=1)
		transcoef=ww[p][k]
		val+=pearson4(transcoef,xx)
	endfor
killwaves transcoef
	return val
end

function TraceMixP4(x1,x2,nbpts, cop4, name)
variable x1,x2,nbpts
string name, cop4
string nameout
make/O/N=(nbpts) transx, transy
transx= min(x1,x2)+abs(x1-x2)*x/nbpts
nameout=name+"_MixP4_outX"
duplicate/O transx $nameout
transy= mixer($cop4,transx[x])
nameout = name+"_MixP4_outY"
duplicate/O transy $nameout
killwaves transx, transy
end

function AddCoP4(int, mass, reso, shape, asym, coefp4)
variable int, mass, reso, shape, asym
string coefp4
make/O/D/N=(5,1) transcoAddCoP4
transcoAddCoP4[][0]={int, mass, reso, shape, asym}
	if(waveExists($coefp4))
		concatenate/NP=1 {transcoAddCoP4}, $coefp4
		killwaves transcoAddCoP4
	else
		rename transcoAddCoP4, $coefp4
	endif
end

function simu2coef(wint,wmass,reso, shape, asym, name)
wave wint, wmass
variable reso, shape, asym
string name
variable k,n=numpnts(wint)
if(waveexists($name))
	killwaves $name
endif
	for(k=0;k<n;k+=1)
		addcop4(wint[k], wmass[k], reso, shape, asym, name)
	endfor
end

function variancep4(reso,shape,asym)
variable reso,shape,asym
variable r=2*(shape-1)
return reso^2*(r^2+asym^2)/(r^2*(r-1))
end

function mixsimu()
wave simu_proba, simu_mass
nvar simu_reso, simu_shape, simu_asym
variable binf,bsup
simu2coef(simu_proba,simu_mass,simu_reso, simu_shape, simu_asym, "simu_cop4")
binf=wavemin(simu_mass)-1
bsup=wavemax(simu_mass)+1
tracemixp4(binf,bsup,500,"simu_cop4","simu_profile")
return sqrt(variancep4(simu_reso, simu_shape, simu_asym))
end


//Revue le 19 février pour ajouter la simulation en profile
function labelize(sima,sipr)//genere les informations masses, abondance relative et deltappm à afficher dans l'elaborateur
wave/Z sima, sipr
wave/Z wave0stripped, deltam
NVAR dismode, simu_prof_sigma
variable maxpro, maxmass , k,n
Duplicate/O sima labelmass
Duplicate/O sipr labelpro
maxpro=wavemax(labelpro)
findvalue/V=(maxpro) labelpro
maxmass=labelmass[V_value]//maxmass est la masse pour laquelle la proba est max
	if (dismode==1)
		labelmass=labelmass-maxmass
	endif
	labelpro=1e6*deltam[x]/sima[x]
simu_prof_sigma=mixsimu()
end

function diffp4(w,x) : fitfunc
wave w
variable x
return pearson4(w,x)* (2*w[1]*w[3]+w[2]*w[4]-2*w[3]*x) / (w[1]^2-2*w[1]*x+w[2]^2+x^2)
end

function inflexp4(reso, shape, asym)
variable reso, shape, asym
return 0.5*reso/shape*sqrt((4*shape^2+asym^2)/(2*shape+1))
end

function setMonoP4(w,np)
wave w
variable np
variable infl=inflexp4(w[2],w[3],w[4]), b=-(infl*(1-10/np))^(1/3), a=((infl*5)^(1/3)-b)
string name
name="Peak_"+num2str(w[1])+"_X"
make/O/D/N=1 transx
transx=w[1]
duplicate/O transx $name
np=2*(floor(np/2))+1
make/D/O/N=(floor(np/2)) transx, transy
transx= w[1]+infl +(a*x/((numpnts(transx)-1))+b)^3
concatenate {transx}, $name
transx= w[1]-infl-(a*x/((numpnts(transx)-1))+b)^3
concatenate {transx}, $name
duplicate/O $name, transx, transy
sort transx, transx, transy
transy= pearson4(w,transx)
duplicate/O transx, $name
name="Peak_"+num2str(w[1])+"_Y"
duplicate/O transy, $name
killwaves transx, transy
end

function textif2cop4(cop4)
string cop4
wave wave1stripped, wave0stripped
variable k, n=numpnts(wave0stripped)
for (k=0;k<n;k+=1)
	AddCoP4(wave1stripped[k], wave0stripped[k], 100, 1000, 0, cop4)
endfor
end

function AddPeakFromAgreg()
getmarquee/W=agregateur bottom
nvar simu_reso,simu_shape, simu_asym
variable n=mass2pnt(V_right)-mass2pnt(V_left), height, center, width=simu_reso, shape=simu_shape, asym=simu_asym
variable tracemin, tracemax
wave roi0
if(stringmatch(S_marqueeWin,"agregateur"))
	Make/O/D/N=(n) roipnts
	roipnts=round(mass2pnt(V_left)+x)
	roipnts2roiz()
	height=wavemax(roi1)
	findvalue/V=(height) roi1
	center=roi0[V_value]
	width=23.6*(V_right-V_left) 
	wave PeakSet_Current
	addcop4(height, center, width, shape,asym, "PeakSet_Current")
	setMonoP4({height, center, width, shape,asym},50)
	addapeaktrace(num2str(center))
	//make/O/D/N=(max(1,dimsize(PeakSet_Current,1))) transpoz, transwidth
	//transpoz=PeakSet_Current[1][x]
	//transwidth=PeakSet_Current[2][x]
	//tracemin=wavemin(transpoz)-0.1
	//tracemax=wavemax(transpoz)+0.1
	//print width/(V_right-V_left)
	//TraceMixP4(tracemin,tracemax,2000, "PeakSet_Current", "PeakResult")
	make/O/D/N=0 roipnts
	roipnts2roiz()
	//killwaves transpoz
endif
SetTraceMixPeaks(PeakSet_Current)
end

function addapeaktrace(cen)
string cen
string tracey, tracex
tracey="Peak_"+cen+"_Y"
tracex="Peak_"+cen+"_X"
Appendtograph/W=agregateur $tracey vs $tracex
ModifyGraph/W=agregateur mode($tracey)=7, marker($tracey)=10,usePlusRGB($tracey)=1;DelayUpdate
ModifyGraph/W=agregateur rgb($tracey)=(0,26112,26112);DelayUpdate
ModifyGraph/W=agregateur plusRGB($tracey)=(0,52224,52224)
ModifyGraph/W=agregateur hbFill($tracey)=30
ReorderTraces/W=agregateur wave1stripped,{$tracey}
end

function removepeaktrace(cen)
string cen
string tracey
tracey="Peak_"+cen+"_Y"
removefromgraph/w=agregateur $tracey
return str2num(cen)
end

function destroypeakwaves(cen)
string cen
string tracey, tracex
tracey="Peak_"+cen+"_Y"
tracex="Peak_"+cen+"_X"
killwaves $tracey, $tracex
return str2num(cen)
end

function majpeaksgraph(ww)//met à jour toutes les waves pour tracer les peaks
wave ww
variable j,k,n=dimsize(ww,1), m=dimsize(ww,0)
for(k=0;k<n;k+=1)
	
endfor
end

function killALLpeaks(ww)
wave ww
variable k,n=max(1,dimsize(ww,1))
make/O/D/N=(n) listcenter
listcenter=ww[1][x]
listcenter=removepeaktrace(num2str(listcenter[x]))
listcenter=destroypeakwaves(num2str(listcenter[x]))
make/O/D/N=(5,0) PeakSet_Current
deletepoints/M=1 0,1, PeakSet_Current
//make/O/N=0 roipnts
//roipnts2roiz()
SetTraceMixPeaks(ww)
killwaves listcenter
end

Menu "graphmarquee"
	"Add this Mass Range to ROI", addROIfromagreg()
	"Substract this Mass Range from ROI", cutROIfromagreg()
	submenu "Peak Management"
	"Add a Peak", addpeakfromagreg()
	"Fit those peaks", FitThosePeaks(PeakSet_Current, hold2cons())
	end
end

function hightlightOnePeak(pozy, pozx, ww)
variable pozy, pozx
wave ww
wave selpklist
variable n, cen
string tracey
string/G hightlightedPeak
variable/g peakpointeur
nvar inco
if(dimsize(ww,0)!=0)
	n=max(1,dimsize(ww,1))
	make/O/D/N=(n) distcenter
	distcenter=abs(pozx-ww[1][x])
	findvalue/V=(wavemin(distcenter)) distcenter
	cen=ww[1][V_value]
	sort distcenter distcenter
	tracey="Peak_"+num2str(cen)+"_Y"
		if(0<pozy && pozy<pearson4({ww[0][V_value],ww[1][V_value],ww[2][V_value],ww[3][V_value],ww[4][V_value]},pozx) && ((n==1)||(distcenter[0]<=distcenter[1]/2)))
			ModifyGraph/W=agregateur rgb($tracey)=(0,26112,0),plusRGB($tracey)=(0,52224,0)
			ModifyGraph/W=agregateur lsize($tracey)=10
			hightlightedPeak=num2str(cen)
			peakpointeur=v_value
			selpklist[v_value][1]=3
		elseif(waveexists($tracey))
			ModifyGraph/W=agregateur rgb($tracey)=(0,26112,26112),plusRGB($tracey)=(0,52224,52224)
			ModifyGraph/W=agregateur lsize($tracey)=1
			hightlightedPeak=""
			peakpointeur=-1
			cen=inco
			selpklist[v_value][1]=2
		endif
else
cen=inco
endif
killwaves distcenter
return cen
end

function KillThisPeak(ww,cen)
wave ww
variable cen
wave/t peaklabels
variable n=max(1,dimsize(ww,1)), m=dimsize(ww,0),  k
findvalue/V=(cen) ww
k=floor(V_value/m)
deletepoints/M=1 k,1,ww
deletepoints k, 1, peaklabels
removepeaktrace(num2str(cen))
destroypeakwaves(num2str(cen))
SetTraceMixPeaks(ww)
return cen
end

function FitThisPeak(ww,cen,cons)
wave ww
variable cen
string cons
variable n=max(1,dimsize(ww,1)), m=dimsize(ww,0),  k
string tracey, tracex
wave wave1stripped
nvar simu_reso,simu_shape, simu_asym
findvalue/V=(cen) ww
k=floor(V_value/m)
tracey="Peak_"+num2str(cen)+"_Y"
tracex="Peak_"+num2str(cen)+"_X"
duplicate/O $tracey cropy
duplicate/O $tracex cropx
cropy=wave1stripped(mass2pnt(cropx))
Make/O/D/N=(m) transco
transco=ww[x][k]
make/O/T/N=2 T_Cons={"K2>0","K3>0"} 
FuncFit/W=2/H=cons/NTHR=0 pearson4 transco cropy /X=cropx /D=$tracey //C=T_Cons
KillThisPeak(ww,cen)
addcop4(transco[0], transco[1], transco[2], transco[3],transco[4], "PeakSet_Current")
setMonoP4(transco,50)
addapeaktrace(num2str(transco[1]))
SetTraceMixPeaks(ww)
simu_reso=transco[2]
simu_shape=transco[3]
simu_asym=transco[4]
killwaves cropy, cropx, transco
end

function FitThosePeaks(ww, cons)//basé sur la pose de marquee, manuellement ou par automatisation
wave ww
string cons
string tracex, tracey,cons2fit
wave wave1stripped
variable n=max(1,dimsize(ww,1)), m=dimsize(ww,0), k
cons2fit=""
getmarquee/W=agregateur bottom
make/o/d/n=(m) transco
make/o/d/n=0 cropx, coef2fit
for(k=0;k<n;k+=1)
	if(ww[1][k]<V_right && ww[1][k]>V_left)
		transco=ww[x][k]
		addcop4(transco[0], transco[1], transco[2], transco[3],transco[4], "trans2fit")
		tracex="Peak_"+num2str(transco[1])+"_X"
		concatenate/NP=0 {$tracex}, cropx
		concatenate/NP=0 {transco}, coef2fit
		cons2fit=cons2fit+cons
	endif
endfor
wave trans2fit
sort cropx cropx
duplicate cropx cropy
cropy=wave1stripped(mass2pnt(cropx[x]))
n=max(1,dimsize(trans2fit,1))
make/o/d/n=(n) cen2fit
cen2fit=trans2fit[1][x]
cen2fit=killthispeak(ww,cen2fit[x])
FuncFit/H=cons2fit/NTHR=0 trans2FitMix coef2fit cropy /X=cropx
for(k=0;k<n;k+=1)
	transco=trans2fit[x][k]
	addcop4(transco[0], transco[1], transco[2], transco[3],transco[4], "PeakSet_Current")
	setMonoP4(transco,50)
	addapeaktrace(num2str(transco[1]))
endfor
SetTraceMixPeaks(ww)
killwaves transco, trans2fit, cen2fit, coef2fit, cropx, cropy
end

function trans2FitMix(w,ix) : fitfunc
wave w
variable ix
wave trans2fit
trans2fit=w
return mixer(trans2fit,ix)
end

function SetTraceMixPeaks(ww)//Met a jour les MixPeaks et SticksPeaks, ne redessine pas les peaks
wave ww
wave/t peaklabels
variable k, n,m
wave wave1stripped, wave0stripped
Make/O/D/N=0 '_MixPeaks_0stripped', '_MixPeaks_1stripped'
string lex, lesx=wavelist("Peak_*_X",";","")
string ley, lesy=wavelist("Peak_*_Y",";","")
n=itemsinlist(lesx)
m=dimsize(ww,0)
make/o/d/n=(n) '_StickPeaks_1stripped', '_StickPeaks_0stripped', '_DataCountsPeaks_1stripped', '_DataCountsPeaks_0stripped'
if(numpnts(ww)!=0)
for(k=0;k<n;k+=1)
	lex=stringfromlist(k,lesx)
	ley=stringfromlist(k,lesy)
	'_StickPeaks_1stripped'[k]=areaXY($lex,$ley)
	'_DataCountsPeaks_1stripped'[k]=areaXY(wave0stripped,wave1stripped,wavemin($lex),wavemax($lex))
	'_StickPeaks_0stripped'[k]=mean($lex)
	'_DataCountsPeaks_0stripped'[k]=mean($lex)
endfor
concatenate/O/NP=0 lesx, '_MixPeaks_0stripped'
sort '_MixPeaks_0stripped' '_MixPeaks_0stripped'
duplicate/O '_MixPeaks_0stripped' '_MixPeaks_1stripped'
'_MixPeaks_1stripped'=mixer(ww,'_MixPeaks_0stripped')
endif
Make/O/T/N=(n,m+2) PeakList
Make/O/T/N=(m+2) wpeaklistlabel
Make/O/T/N=(n) peaklabels
PeakList=num2str(ww[y][x])
peaklist[][m]=peaklabels[x]
peaklist[][m+1]=""
Make/O/N=(n,m+2) SelPkList=2
selpklist[][m+1]=32
SetResidueMixPeaks(ww)
listelesdonnees()
end

function SetResidueMixPeaks(ww)
wave ww
wave '_MixPeaks_0stripped', '_MixPeaks_1stripped', wave1stripped
duplicate/O '_MixPeaks_1stripped' MixPeaksResidue_Y
MixPeaksResidue_Y= 100*('_MixPeaks_1stripped' - wave1stripped[mass2pnt('_MixPeaks_0stripped')])/wavemax('_MixPeaks_1stripped')
if (numpnts(ww)==0)
	Make/o/d/n=0 listzentrum, listwidz
else
	matrixop/O listzentrum=(row(ww,1))^t
	matrixop/O listwidz=(row(ww,2))^t
endif
make/o/T/n=(numpnts(listzentrum)) labelpeaks
labelpeaks=num2str(listzentrum[x])
end

Window PeakQualCont() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(754.5,42.5,1335.75,461)/K=1 MixPeaksResidue_Y vs '_MixPeaks_0stripped'
	Movewindow/W=PeakQualCont 830.25, 42.5, 1389, 484.25
	AppendToGraph/T MixPeaksResidue_Y vs '_MixPeaks_0stripped'
	AppendToGraph/L=axemix '_MixPeaks_1stripped' vs '_MixPeaks_0stripped'
	AppendToGraph/R=axestix '_StickPeaks_1stripped' vs '_StickPeaks_0stripped'
	ModifyGraph userticks(top)={listzentrum,labelpeaks}
	ModifyGraph margin(right)=227
	ModifyGraph mode(MixPeaksResidue_Y)=7,mode(MixPeaksResidue_Y#1)=3,mode('_StickPeaks_1stripped')=1
	ModifyGraph marker(MixPeaksResidue_Y#1)=62
	ModifyGraph lSize('_StickPeaks_1stripped')=2
	ModifyGraph rgb(MixPeaksResidue_Y)=(0,0,0),rgb(MixPeaksResidue_Y#1)=(0,0,0),rgb('_MixPeaks_1stripped')=(0,26112,26112)
	ModifyGraph rgb('_StickPeaks_1stripped')=(0,52224,52224)
	ModifyGraph msize(MixPeaksResidue_Y#1)=1
	ModifyGraph hbFill(MixPeaksResidue_Y)=31
	ModifyGraph usePlusRGB(MixPeaksResidue_Y)=1
	ModifyGraph useNegRGB(MixPeaksResidue_Y)=1
	ModifyGraph hBarNegFill(MixPeaksResidue_Y)=31
	ModifyGraph plusRGB(MixPeaksResidue_Y)=(0,52224,0)
	ModifyGraph grid(top)=1
	ModifyGraph noLabel(axestix)=2
	ModifyGraph fStyle(left)=1,fStyle(bottom)=1,fStyle(axemix)=1,fStyle(axestix)=1
	ModifyGraph axThick(left)=2,axThick(bottom)=2,axThick(top)=2,axThick(axemix)=2,axThick(axestix)=0
	ModifyGraph gridRGB(top)=(47872,47872,47872)
	ModifyGraph gridStyle(top)=5
	ModifyGraph lblPos(left)=67
	ModifyGraph tkLblRot(top)=90
	ModifyGraph freePos(axemix)=0
	ModifyGraph freePos(axestix)=0
	ModifyGraph axisEnab(left)={0,0.45}
	ModifyGraph axisEnab(axemix)={0.55,1}
	ModifyGraph axisEnab(axestix)={0.55,1}
	Label left "Not fit intensity rel. to max. peak (%)"
	Label bottom "m/z"
	SetAxis left -100,100
	SetAxis bottom 96.6065796440493,98.8365705649496
	SetAxis top 96.6065796440493,98.8365705649496
	SetAxis axemix 0,*
	SetAxis axestix 0,*
	ControlBar 114
	ListBox list0,pos={10,22},size={84,85},fSize=6,listWave=root:wpeaktohold
	ListBox list0,selWave=root:sel_peaktohold,mode= 5
	GroupBox group0,pos={3,3},size={98,109},title="Hold:"
	ListBox list1,pos={465,122},size={297,715},proc=ListBoxPeakControl
	ListBox list1,listWave=root:PeakList,selWave=root:SelPkList
	ListBox list1,titleWave=root:wpeaklistlabel,mode= 5,editStyle= 1
	ListBox list1,userColumnResize= 1
	Button button0,pos={648,5},size={107,20},proc=killtickedpeaks_butt,title="Delete ticked peaks"
	Button button0,fColor=(65280,0,0)
	Button producebaseline,pos={109,4},size={202,20},proc=PruduceBaselineButton,title="Produce baseline from peaks feet"
	Button producebaseline,fColor=(65280,43520,0)
	Button button1,pos={648,26},size={107,20},proc=tickallbutton,title="Tick all"
	Button button1,fColor=(65280,43520,0)
	Button button2,pos={648,48},size={107,20},proc=fittickedbutton,title="Fit ticked"
	Button button2,fColor=(65280,0,0)
EndMacro

function/S hold2cons()
wave sel_peaktohold
variable k, n=numpnts(sel_peaktohold)
string res=""
for(k=0;k<n;k+=1)
	res+=num2str(sel_peaktohold[k]==48)
endfor
return res
end

function calculatefrompeak()
nvar inco, findindata, peakpointeur
svar nommol
wave/Z list_molperm, molecule, incoy
wave/T peaklabels
string lawave
variable transfindindata=findindata
		if(numpnts(list_molperm)>0)
			analysemass("arg")
			setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
			setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
			lawave=stringfromlist(0,wavelist("molprop_*",";",""))
			molecule=0
			tuelesspec()
			duplicate/O $lawave molecule
			findindata=2
			genestringmol()
			peaklabels[peakpointeur]=nommol
			findindata=transfindindata
			ListBox list0 win=attributeur, selRow=0
			SetTraceMixPeaks(PeakSet_Current)
		else
			razprop()
		endif
end

function AddPeakFromCurrentSimu()
wave simu_proba, simu_mass
nvar simu_reso, simu_shape, simu_asym
variable k, n=numpnts(simu_proba)
for(k=0;k<n;k+=1)
	addcop4(simu_proba[k], simu_mass[k], simu_reso, simu_shape, simu_asym, "PeakSet_Current")
	setMonoP4({simu_proba[k], simu_mass[k], simu_reso, simu_shape, simu_asym},50)
	addapeaktrace(num2str(simu_mass[k]))
endfor
SetTraceMixPeaks(PeakSet_Current)
end

Function ListBoxPeakControl(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
wave PeakSet_Current
variable n,m
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
m=dimsize(selwave,1)-2
n=max(1,dimsize(selwave,0))
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 3: // double click
			break
		case 4: // cell selection
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			make/O/D/N=(m) transco
			transco=str2num(listwave[row][x])
			//print PeakSet_Current[1][row]
			ReplacePeak(PeakSet_Current,transco, PeakSet_Current[1][row])
			killwaves transco
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch

	return 0
End

function FindColFromCen(ww,cen)
variable cen
wave ww
variable n=max(1,dimsize(ww,1)),m=dimsize(ww,0)
switch(numpnts(ww))
	case 0:
		return 0
	default:
		findvalue/V=(cen) ww
		return floor(V_value/m)
endswitch
end

function ReplacePeak(ww,w, oldcen)
wave ww, w
variable oldcen
wave/T peaklist, peaklabels
variable k
string oldnamex, newnamex, oldnamey, newnamey
switch(numpnts(ww))
	case 0:
		return 0
	default:
		k=FindColFromCen(ww,oldcen)
		ww[][k]=w[x]
		oldnamex="Peak_"+num2str(oldcen)+"_X"
		newnamex="Peak_"+num2str(w[1])+"_X"
		oldnamey="Peak_"+num2str(oldcen)+"_Y"
		newnamey="Peak_"+num2str(w[1])+"_Y"
		rename $oldnamey, $newnamey
		rename $oldnamex, $newnamex
		setMonoP4(w,50)
		peaklabels[k]=peaklist[k][5]
		SetTraceMixPeaks(ww)
endswitch
end

/////

//ajout aux fonctions de gestion des molbags
function MolbagKill(nom)
string nom
string nommolbag, nomchargebag, nomlistbag, nomciblebag
nommolbag="Molbag_"+nom
nomchargebag="Chargebag_"+nom
nomlistbag="Listbag_"+nom
nomciblebag="Ciblebag_"+nom
killwaves/Z $nommolbag, $nomchargebag, $nomlistbag, $nomciblebag
end

function perm2current()
string lesperm=wavelist("molperm_*",";",""), laperm
nvar charge
nvar krikriter
variable k, n=itemsinlist(lesperm), transcharge=charge, transkrikriter=krikriter
wave molecule
charge=0
krikriter=1
for(k=0;k<n;k+=1)
	laperm=stringfromlist(k,lesperm)
	molecule=0
	tuelesspec()
	duplicate/O $laperm molecule
	genestringmol()
	addthissimu()
endfor
charge=transcharge
krikriter=transkrikriter
end

function current2perm()
string lesmol=wavelist("Molecule_*",";",""), lamol
variable k, n=itemsinlist(lesmol)
wave molecule
for(k=0;k<n;k+=1)
	lamol=stringfromlist(k,lesmol)
	molecule=0
	tuelesspec()
	duplicate/O $lamol molecule
	genestringmol()
	addmolTOperm()
endfor
end

function emptyPerm()
string lesmol=wavelist("Molperm_*",";",""), lamol
variable k, n=itemsinlist(lesmol)
for(k=1;k<(n+1);k+=1)
	killthismol(0)
endfor
end

Function ButtonemptyPerm(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			emptyPerm()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Buttoncurrent2perm(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			current2perm()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			addthemolbag(popstr)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



Function PopMenuSavetoMolbag(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
string name
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			strswitch(popstr)
				case "_New":
					prompt name, "Name for new molecule set"
					doprompt "Name for new molecule set", name
					fullsavebag(name)
					break
				case "Factors":
					fullsavebag(popstr)
					emptyperm()
					current2perm()
					calcperm()
					break
				default :
					fullsavebag(popstr)
					break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function PopMenuKillMolbag(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
string yorn
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			prompt yorn, "Are you sure ?", popup, "Yes;No"
			Doprompt "Are you sure ?", yorn
			strswitch(yorn)
			case "yes":
			MolbagKill(popstr)
			break
			case "No" :
			break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function transformcurrentset(multiplier)
variable multiplier
wave molecule
nvar charge
fullsavebag("transbag")
duplicate/O molecule transmol
wave Molbag_transbag, chargebag_transbag
variable k, n=max(1,dimsize(Molbag_transbag,1))
for(k=0;k<n;k+=1)
	Molbag_transbag[][k]+=transmol(x)*multiplier
endfor
chargebag_transbag=charge
videagreg()
addthemolbag("transbag")
molecule=0
tuelesspec()
duplicate/O transmol molecule
genestringmol()
MolbagKill("transbag")
killwaves transmol
end

Function MolbagOperateToALL(bag,molecule,charge)
string bag
wave molecule
variable charge
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1),k
chargebag=charge
for(k=0;k<n;k+=1)
	molbag[][k]+=molecule[x]
	extract/O/FREE molbag,transmol, y==k
	listbag[k]=mol2str(transmol,charge)
	ciblebag[k][1]=mol2massmono(transmol,charge)
	ciblebag[k][3]=0
endfor
end

Function ButtonSumCurrentFormula(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			transformcurrentset(1)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonSubCurrentFormula(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			transformcurrentset(-1)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//////tri sur segment, autocorrelation isotopes, statistiques sur intercepts
function autocorr(centre, gamme, pas)
variable centre, gamme, pas
wave wave1stripped, wave0stripped
variable t, k,n=numpnts(wave1stripped)
string namey, namex
namey="AutocorrY_"+num2str(centre)
namex="AutocorrX_"+num2str(centre)
Make/O/D/N=(2*gamme/pas) transy, transx
transx=(centre-gamme)+x*pas
transy=0
for(k=0;k<n;k+=1)
	transy+=wave1stripped[k]*wave1stripped[mass2pnt(wave0stripped[k]+transx[x])]
	//doupdate
endfor
//Duplicate/O transy, $namey
//Duplicate/O transx, $namex
//killwaves transy, transx
end

function tri2lincoord(row,col,n,k)
variable row,col,n,k
//variable k=edge2numpnts(n)
return (k-edge2numpnts(n-1-col)-(n-row))*(row>col)
end

function x2tricol(p,k)
variable p,k
variable n=numpnts2edge(k), q=k-p
return n-ceil(numpnts2edge(q))
end

function x2trirow(p,k)
variable p,k
variable q=k-p, n=numpnts2edge(k)
return n-(q-edge2numpnts(floor(numpnts2edge(q-1))))
end

function numpnts2edge(k)//donne la dimension du carré qui contient une triangulaire stricte de k points
variable k
return (sqrt(8*k+1)+1)/2
end

function edge2numpnts(n)//donne le nombre de point de la triangulaire stricte issue de n*n
variable n
return (n*n-n)/2
end

function filltrig(ww)
wave ww
variable n=dimsize(ww,0), m= dimsize(ww,1), i ,j, k=0
if(n==m)
	for(j=0;j<n;j+=1)
		for(i=j+1;i<n;i+=1)
			ww[i][j]=k
			k+=1
		endfor
	endfor
endif
end

function compareboucle(n)
variable n
make/O/D/N=(n,n) testmat, testmat2
testmat=0
testmat2=0
variable t1,t2, k =edge2numpnts(n)
t1=ticks
filltrig(testmat)
t1=ticks-t1
t2=ticks
testmat2=tri2lincoord(x,y,n,k)
t2=ticks-t2
print t1/60,t2/60
end

/////////////////////////////////////////////////////////
/// MAJ Roland Thissen 23 avril 2014 ///
/////////////////////////////////////////////////////////

function get_generalHook(info)
string info
string lawave, lordre, lepoint
string molecule
wave element_H
wave element_C
wave element_N
wave element_O

if(stringmatch(stringbykey("event",info),"mouseup")==1)
make/o/d/n=1 pozx, pozy, controlenfonce
pozx=str2num(stringbykey("mousex",info))
pozy=str2num(stringbykey("mousey",info))
controlenfonce=str2num(stringbykey("MODIFIERS",info))

	if(strlen(TraceFromPixel (pozx[0],pozy[0],""))!=0)							// 18 avril 2014 RT elimine le bug quand on tape trop loin d'un point
		lepoint=stringbykey("HITPOINT",TraceFromPixel (pozx[0],pozy[0],"")) 
		dowindow/F vankrevH_N
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,65280,0) A "+ lawave+ " " + lepoint
		execute lordre 
		dowindow/F GrN_C
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,0,0) A "+ lawave+ " " + lepoint
		execute lordre  
		dowindow/F GrH_C
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,0,0) A "+ lawave+ " " + lepoint
		execute lordre 
		dowindow/F GrDBE
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,0,0) A "+ lawave+ " " + lepoint
		execute lordre 
		molecule="\\Z14\\f01C\\B"+num2str(element_C[str2num(lepoint)])
		
		molecule+="\\M\\Z14H\\B"+num2str(element_H[str2num(lepoint)])
		if (element_N[str2num(lepoint)]>0)
			molecule+="\\M\\Z14N\\B"+num2str(element_N[str2num(lepoint)])
		endif
		if (element_O[str2num(lepoint)]>0)
			molecule+="\\M\\Z14O\\B"+num2str(element_O[str2num(lepoint)])
		endif
		TextBox/C/N=text0/F=0/A=MC/X=-44.21/Y=45.33 molecule

		dowindow/F tholinomics
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,65280,0) A "+ lawave+ " " + lepoint
		execute lordre 
		dowindow/F Grcomposite
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,65280,0) A "+ lawave+ " " + lepoint
		execute lordre 
		dowindow/F GrZubarev
		lawave=WaveName("",0,1)
		lordre="cursor/S=2/C=(0,0,0) A "+ lawave+ " " + lepoint
		execute lordre 
		if (wavemax(element_o)>0)
			dowindow/F vankrevO
			lawave=WaveName("",0,1)
			lordre="cursor/S=2/C=(0,65280,0) A "+ lawave+ " " + lepoint
			execute lordre 
		endif
		if (controlenfonce(0)>=8)
			removecursorpoint("bob")
		endif
	endif
endif
return 1
end

////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
///////////////// UPGRADES FOR GLURNS ///////
///////////////// 24 avril 2014 /////////////////////////////

Function choisyladata(lba) : ListBoxControl //////// OBSOLETE
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listwave= lba.listWave
	WAVE/Z selwave = lba.selWave

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 3: // double click
			break
		case 4: // cell selection
		activeladata("")
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
	endswitch

	return 0
End

Function ActiveLaDataButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			controlinfo/W=advancedmanager TabCon_List_DATA
			if(v_value==-1)
				break
			endif
			wave/T listwave=$S_value
			string name
			name=listwave[V_value][1]
			activeladata(name)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


function LoadCosima()// !!! cosimadata0 est l'intnesité
    string cheminversfichier
    string nomfichier, fichiery, fichierx, fichierinfo
    GetFileFolderInfo/Q
    NewPath/O/Q ledossier ParseFilePath(1, S_Path, ":", 1, 0)
    nomfichier=ParseFilePath(0, S_Path, ":", 1, 0)
    LoadWave/A=CosimaHeader/O/J/K=2/Q/D/L={0,0,93,0,0}/P=ledossier nomfichier   
    LoadWave/A=CosimaData/O/G/Q/D/L={0,0,0,1,2}/P=ledossier nomfichier   
    print nomfichier
    wave/T cosimaheader0
    wave cosimadata0, cosimadata1
    cosimaheader0=replacestring(" ",cosimaheader0,"")
    cosimaheader0=replacestring("\"",cosimaheader0,"")
    string target="",coordx="",coordy="", name="", prodid=""
    variable n=numpnts(cosimaheader0),k
    for(k=0;k<n;k+=1)
        target+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_ID",cosimaheader0[k], "=")
        coordx+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_X",cosimaheader0[k], "=")
        coordy+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_Y",cosimaheader0[k], "=")
        prodid+=stringbykey( "PRODUCT_ID",cosimaheader0[k], "=")
    endfor
    name=prodid+"_"+target+"_"+coordx+"_"+coordy
    name=name[3,22]
    string loadoption
    prompt name, "Name this Sample"
    prompt loadoption, "On-the-fly loading options", popup, "Load it as it is;Delete 0s intensities;Keep >0.5 masses"
    doprompt "Enter a good name for this sample",name, loadoption
    if(v_flag==0)
        strswitch(loadoption)
            case "Load it as it is":
               
                break
            case "Delete 0s intensities":
                duplicate/O cosimadata0 cropcosi
                cropcosi=cosimadata0/cosimadata0*x
                wavetransform zapnans cropcosi
                duplicate/O cropcosi, cropcosi0, cropcosi1
                cropcosi0=cosimadata0[cropcosi[x]]
                cropcosi1=cosimadata1[cropcosi[x]]
                duplicate/O cropcosi0 cosimadata0
                duplicate/O cropcosi1 cosimadata1
                killwaves cropcosi, cropcosi0, cropcosi1
                nomfichier+="z"
                break
            case "Crop only >1 intensities":
                break
            case "Keep >0.5 masses" :
                deletepoints 0,1, cosimadata0, cosimadata1   
                findlevel/Q cosimadata1, 0.5
                deletepoints 0,V_levelx, cosimadata0, cosimadata1
                cosimadata0=max(cosimadata0,0)
                break
        endswitch
    nomfichier=name[0,19]
    fichierinfo=nomfichier+"_info"
    fichiery=nomfichier+"_1stripped"
    fichierx=nomfichier+"_0stripped"
    Rename cosimadata1, $fichierx
    Rename cosimadata0, $fichiery
    Rename cosimaheader0, $fichierinfo
    else
    killwaves cosimaheader0, cosimadata1, cosimadata0
    endif
    listelesdonnees()
end

function AddPeakFromROI()
nvar simu_reso,simu_shape, simu_asym
variable height, center, width=simu_reso, shape=simu_shape, asym=simu_asym
wave roi0, roi1, roi2
wave PeakSet_Current
variable n=numpnts(roi0),k
for(k=0;k<n;k+=1)
height=roi1[k]
center=roi0[k]
addcop4(height, center, width, shape,asym, "PeakSet_Current")
setMonoP4({height, center, width, shape,asym},50)
addapeaktrace(num2str(center))
endfor
SetTraceMixPeaks(PeakSet_Current)
end

function killtickedpeaks()
wave selpklist, peakset_current
wave/t peaklist
variable n= dimsize(selpklist,0),k
for(k=n-1;k>-1;k-=1)
	if(selpklist[k][6]==48)
		print k, peaklist[k][1]
		deletepoints/M=1 k,1,peakset_current
		deletepoints k, 1, peaklabels
		removepeaktrace(peaklist[k][1])
		destroypeakwaves(peaklist[k][1])
	endif
endfor
SetTraceMixPeaks(PeakSet_Current)
end

Function killtickedpeaks_butt(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			killtickedpeaks()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function switch_profile(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
nvar simu_shape, simu_reso
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			switch(popnum)	// numeric switch
				case 1:
					simu_reso=10
					simu_shape=1e6
					break
				case 2:
					simu_reso=0.01
					simu_shape=1
					break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function bouton1(NomButton) : ButtonControl
	String NomButton
	variable ufparam, radioparam
	controlinfo/W=filtres setvar0
	ufparam=v_value
	controlinfo/W=filtres setvar1
	radioparam=v_value
	detectionUF_et_filtreRadio(wave0stripped, wave1stripped,ufparam, radioparam, "ROI")	
	calcinfo()
End

Function updateUF(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			variable ufparam, radioparam
			controlinfo/W=filtres setvar0
			ufparam=v_value
			controlinfo/W=filtres setvar1
			radioparam=v_value
			detectionUF_et_filtreRadio(wave0stripped, wave1stripped,ufparam, radioparam, "ROI")
			calcinfo()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function UpdateRadio(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			variable ufparam, radioparam
			controlinfo/W=filtres setvar0
			ufparam=v_value
			controlinfo/W=filtres setvar1
			radioparam=v_value
			detectionUF_et_filtreRadio(wave0stripped, wave1stripped,ufparam, radioparam, "ROI")
			calcinfo()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function ButtonProc_5(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=advancedmanager TabCon_List_DATA
			wave/T listwave=$S_value
			string newname=listwave[v_value][1]+"_ROIed"
			prompt newname, "New Name"
			Doprompt "Enter a new good name for this sample",newname
			newname=newname[0,19]
			if(v_flag==0)
				roi2data(newname)
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


///////////Baseline stripper cubic cpline

/// remettre la liste des points à 0
function razbaseline()
make/O/N=0 basex, basey
end

/// maintenir la liste des points et générer la baseline
function ProdBaseline()
wave wave1stripped, basey, basex
duplicate/O wave1stripped, '_baseline_1stripped', '_baseline_0stripped'
Interpolate2/T=3/N=200/E=1/I=3/Y='_baseline_1stripped'/X='_baseline_0stripped' basex, basey
listelesdonnees()
end

/// générer des points à n sigma des pics en évitant les zones contigues...
function localbgrd2baseline(ww,average)
variable average //nb de points pour la moyenne au piede de chaque lot de pics
wave ww
wave basex, basey, '_MixPeaks_0stripped', wave1stripped, wave0stripped
variable n,k
Differentiate '_MixPeaks_0stripped'/D='_MixPeaks_0stripped_DIF'
findlevels/Q/D=localbgrdx '_MixPeaks_0stripped_DIF' 0.2 // ATTENTION paramètre arbitraire de distance entre les points dits contigus
killwaves '_MixPeaks_0stripped_DIF'
insertpoints 0,2, localbgrdx
localbgrdx[0]=0
localbgrdx[1]=numpnts('_MixPeaks_0stripped')-1
localbgrdx='_MixPeaks_0stripped'[localbgrdx]
sort localbgrdx localbgrdx
duplicate/O localbgrdx localbgrdy
n=numpnts(localbgrdx)
for(k=0;k<n;k+=1)
	findlevel/Q wave0stripped localbgrdx[k]
	localbgrdx[k]=mean(wave0stripped, v_levelx-average, v_levelx+average)
	localbgrdy[k]=mean(wave1stripped, v_levelx-average, v_levelx+average)
endfor
concatenate/NP {localbgrdx}, basex
concatenate/NP {localbgrdy}, basey
end

Function PruduceBaselineButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			razbaseline();localbgrd2baseline(peakset_current,1);ProdBaseline()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonKillFatNoise(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave/Z wave1stripped,wave0stripped, roi2, roi0
	switch( ba.eventCode )
		case 2: // mouse up
			deletenoiseFAT(wave1stripped,wave0stripped)
			duplicate/O roi0 roi2
			roi2=defmass(roi0)
			calcinfo()
			break
	endswitch

	return 0
End

Function ButtonKillFineNoise(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave/Z wave1stripped,wave0stripped, roi2, roi0
	switch( ba.eventCode )
		case 2: // mouse up
			deletenoiseFine(wave1stripped,wave0stripped)
			duplicate/O roi0 roi2
			roi2=defmass(roi0)
			calcinfo()
			break
	endswitch

	return 0
End


function autofit()
localbgrd2baseline(peakset_current,1)
end

function copydata2roi()// ATTENTION 100 PARAMETRE ARBITRAIRE
nvar inco, intenscumthres
controlinfo/W=advancedmanager TabCon_List_DATA
wave/T listwave=$S_value
if(v_value==-1)
	return -1
endif
string laun, lazero, latrois
wave/Z topbrighty, topbrightx, incoy, incox, logwave1, incoxdm
laun=listwave[V_value][1]+"_1stripped"
lazero=listwave[V_value][1]+"_0stripped"
latrois=listwave[V_value][1]+"_3stripped"
Duplicate/O $laun roi1
Duplicate/O $lazero roi0
Duplicate/O $lazero roi2
if(waveexists($latrois))
	Duplicate/O $latrois roi3
endif
roi2=defmass(roi0)
Sort roi0, roi0, roi1, roi2, roi3
calcinfo()
end 

Function loadInRoi(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			copydata2roi()
		case -1: // control being killed
			break
	endswitch

	return 0
End

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Djevahirdjian Léo
// été 2014
//
//		Compare Attrib
//			Interface de comparaison et d'étude des résultats d'une
//	attribution, dans une interface qui se veut adaptative basée sur
//	l'utilisation des MolBags.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Compare Attrib
//
//		I. Parties Obligatoires a intégrer dans Attributor
////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compare Attrib
//
//	I.a Anciennes parties d'Attributor à modifier
//		Modification du panneau d'attribution pour activer ou non la
//	visualisation permettant un gain de temps IMMENSE lors des
//	attributions. Remodelages de la fonction d'attribution pour la
//	nouvelle interface d'étude des résultats de celle-ci, plus modulable
//	que celle précédente.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// panel3()
//	Nouveau panneau d'attribution permetant l'affichage ou non des
//	attributions concluantes sur le spectre et les diagrammes dm vs m.
//	Desactiver cetrte option permet un gain considérable, en effet on a
//	une dependance de l'affichage qui semble être en N2 du nombre de
//	graph a afficher dans une même fenêtre.
////////////////////////////////////////////////////////////////////////////////////////////////////////

Window panel3() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /K=1 /W=(188,425,418,595)/FLT
	SetVariable setvar0,pos={3,3},size={224,16},title="Next m/z positions to analyse :"
	SetVariable setvar0,limits={1,inf,1},value= nexttoatt,live= 1
	Button button0,pos={3,126},size={50,40},proc=ButtonGoForAttribution,title="Go !"
	Button button0,fColor=(65280,0,0)
	SetVariable setvar1,pos={27,22},size={200,16},title="For each, check these first"
	SetVariable setvar1,value= nbprop,live= 1
	SetVariable setvar2,pos={27,42},size={200,16},title="With this threshold (ppm)"
	SetVariable setvar2,value= ppmthreshold,live= 1
	SetVariable setvar3,pos={27,62},size={200,16},title="With this ion charge"
	SetVariable setvar3,value= charge,live= 1
	SetVariable setvar4,pos={27,82},size={200,16},title="And this combo probability"
	SetVariable setvar4,limits={0,1,0.01},value= krikriter,live= 1
	CheckBox checkAffichage,pos={5,106},size={199,14},proc=pilotageAffichageAttrib,title="Affichage pendant Attribution ON/OFF"
	CheckBox checkAffichage,variable= affichageAttribONOFF
	Button button1,pos={54,126},size={100,40},proc=GoForAttribROI,title="Go using ROI\ras guide"
	Button button1,fColor=(65280,0,0)
	Button button2,pos={154,126},size={73,40},proc=ButtonFasttribution,title="Fasttribution\r of the ROI"
	Button button2,fColor=(65280,43520,0)
EndMacro

////////////////////////////////////////////////////////////////////////////////////////////////////////
// pilotageAffichageAttrib(nomCheckBox,active)
//	Control de la check box permetant de desactiver ou non l'affichage
//	necessaire à panel3, pilotte une variable globale qui sera lue au
//	moment de l'attribution
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function pilotageAffichageAttrib(nomCheckBox,active) : CheckBoxControl
	String nomCheckBox
	Variable active
	Variable /G affichageAttribONOFF = active
End //  pilotageAffichageAttrib(nomCheckBox,active)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// addthissimu()
// 	Fonction qui gère l'ajout de la molécule attribuée dans les différentes
//	instances necessaires (waves de travail, agregateur, affichage dans
//	spectre et defaut masse vs masse.
//	C'est dans cette fonction que l'affichage a lieu, elle a donc été modifiée.
//	Une condition de pilotage a été ajoutée pour allumer ou éteindre
//	l'affichage des attributions concluantes.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function addthissimu()
	SVAR nommol
	NVAR affichageAttribONOFF, lastincox, lastincoy,charge

	wave/Z simu_proba, simu_mass, Index_formules, wave1stripped, wave0stripped, molecule, List_targets, deltam, test1
	wave/T list_formules
	string blazmass, blazproba,labname="prout", blazmol, blazdefmass
	variable nb_mol=itemsinlist(wavelist("Formule_*",";",""))/2+1, pmaxproba
	variable res=(wavemax(wave0stripped)-wavemin(wave0stripped))/numpnts(wave0stripped)
	Duplicate/O test1, testcrit

	findvalue/V=(wavemax(simu_proba)) simu_proba
	pmaxproba=V_value
	InsertPoints numpnts(Index_formules),1, Index_formules
	Index_formules[numpnts(Index_formules)-1]= simu_mass[pmaxproba]
	InsertPoints numpnts(list_formules),1, List_formules, List_targets
	list_formules[numpnts(list_formules)-1]=nommol
	List_targets[numpnts(list_formules)-1][1]=lastincox
	List_targets[numpnts(list_formules)-1][2]=lastincoy
	List_targets[numpnts(list_formules)-1][3]=1e6*(simu_mass[pmaxproba]-lastincox)/simu_mass[pmaxproba]
	testcrit=test1*simu_proba* abs( List_targets[numpnts(list_formules)-1][3]- 1e6*deltam/simu_mass)
	List_targets[numpnts(list_formules)-1][0]=mean(testcrit)
	
	MolbagAddTo(molecule,charge,nommol,{List_targets[numpnts(list_formules)-1][0],List_targets[numpnts(list_formules)-1][1],List_targets[numpnts(list_formules)-1][2],List_targets[numpnts(list_formules)-1][3],List_targets[numpnts(list_formules)-1][4]},"Current")
	
	variable rouge=0, bleu=3000*variance(molecule),vert=1000
	blazmol= "Molecule_"+num2str(nb_mol)
	blazmass= "Formule_"+num2str(nb_mol)+"mass"
	blazproba="Formule_"+num2str(nb_mol)
	blazdefmass="DefmassForm_"+num2str(nb_mol)
	Duplicate/O molecule $blazmol
	Duplicate/O simu_proba $blazproba
	Duplicate/O simu_mass $blazmass
	Duplicate/O simu_mass simu_defmass
	simu_defmass=simu_mass-round(simu_mass)
	Duplicate/O simu_defmass $blazdefmass
	If(affichageAttribONOFF)
		AppendToGraph/W=agregateur /C=(rouge,bleu,vert), $blazproba vs $blazmass
		ModifyGraph/W=agregateur mode($blazproba)=8, marker($blazproba)=18,useMrkStrokeRGB($blazproba)=1
		AppendToGraph/W=dmvm /C=(rouge,bleu,vert), $blazdefmass vs $blazmass
		ModifyGraph/W=dmvm mode($blazdefmass)=3, marker($blazdefmass)=18, useMrkStrokeRGB($blazdefmass)=1
	endIf
end

////////////////////////////////////////////////////////////////////////////////////////////////////////
// removeFakes()
//	Fonction appelée après l'attribution pour faire un peu de ménage.
//	Elle a été réécrite pour s'adapter à la nouvelle idée de compare
//	attrib plus flexible.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function removeFakes()
	Wave Proximite
	Wave Massinconnue
	Wave Delta
	Wave intensitePic
	
	Wave /T compareEnCoursElementsWaves
	Wave /T listbagCompareEnCours
	
	Duplicate /O Proximite OLDProximite
	Sort OLDProximite, MassInconnue, Delta, intensitePic, Proximite,listbagCompareEnCours
	Variable o
	For(o=0 ; o<numpnts(compareEnCoursElementsWaves) ;o+=1)
		Sort OLDProximite, $compareEnCoursElementsWaves[o]
	endFor
      
      variable niveau
      Findlevel/Q proximite 1
      if (v_flag==0)
      		niveau=v_levelx
             Make/O/N=(numpnts(proximite)-niveau), massegiclee, intensitegiclee
            	massegiclee= massinconnue (niveau+x)
             intensitegiclee=intensitepic (niveau+x)
             Variable niveauEnd = numpnts(Proximite)
             DeletePoints niveau, niveauEnd, MassInconnue, Delta, intensitePic, Proximite,listbagCompareEnCours
             Variable m
		For(m=0 ; m<numpnts(compareEnCoursElementsWaves) ;m+=1)
			DeletePoints niveau, niveauEnd, $compareEnCoursElementsWaves[m]
		endFor
  	endif
      sort massegiclee, massegiclee, intensitegiclee
      goods2data("bel");gicle2data("gic");  
      
      KillWaves OLDProximite
end

////////////////////////////////////////////////////////////////////////////////////////////////////////
// fullauto(nprochains,seuil)
//	Fonction principale, effectuant l'attribution. Elle place le résultat
//	dans le molbag Current et fichiers associés. Modification pour
//	prendre en compte le nouveau format de molbag, plus complet.
//	Execute en fin l'appel pour l'interface compare Attrib.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function fullauto(nprochains,seuil)// ATTENTION POSSIBLES PROBLEMES DU AU FAIT QUE NBPROP PEUT ETRE SUPERIEUR A J, CONTOURNEMENT PAR STR2NUM !!!

	Variable nprochains,seuil
	Variable k, n=nprochains, j

	Nvar inco, nbprop
	
	videagreg()

	// recup des waves du spectre de masse
	Wave/Z wave0stripped, wave1stripped
	// Wave de l'inconu de la mort qui tue
	Wave/Z dataindex,deltam, deltai, molecule, transindex, simu_proba, simu_mass, deltappm, test1, labelpro, incox,incoxdm,incoy, testlocal
	//Waves Text
	Wave/T list_prop, elements,List_formules, settings
	String lawave


	// SI il existe on détruit le panel3
	If(Wintype("panel3")!=0)
		Killwindow panel3
	endIf

	// La wave Progressbarre va contenir toutes les infos
	Make/O/N=(n,118)/D progressbarre
	progressbarre=1e12
	//execute "showprogbar()"
		for(k=0;k<n;k+=1)
			analysemass("arg")
			Make/O/D/N=(nbprop) critere
			critere=1
				for(j=0;j<nbprop;j+=1)
					molecule=0
					tuelesspec()
					lawave="molprop_"+num2str(j)
					duplicate/O $lawave Molecule
					genestringmol()
					duplicate/O test1, testlocal
					testlocal= test1 * simu_proba[x] * abs( str2num(list_prop[j][1])- 1e6*deltam[x]/simu_mass[x])
					duplicate/O simu_proba transindex
					makeindex/R simu_proba transindex
						if(abs(labelpro[transindex[0]]-str2num(list_prop[j][1]))<0.1 && abs(str2num(list_prop[j][1]))<seuil && abs(labelpro[transindex[0]])<seuil)//0.1 arbitraire
							critere[j]=mean(testlocal)
						else
							critere[j]=1
						endif
				endfor
			// recup Proximite
			progressbarre[k][0]= wavemin(critere)
			// recup Masse inconnue
			progressbarre[k][1]= inco
			// recup de l'intensite de la masse inconnue
			progressbarre[k][2]= incoy[0]
			if(wavemin(critere)<1)
				Findvalue/V=(wavemin(critere)) critere
				molecule=0
				tuelesspec()
				lawave="molprop_"+num2str(v_value)
				duplicate/O $lawave Molecule
				// genere nom molecule
				genestringmol()
				//ajoutte affichage de masse trouvée et ajoute dans liste pour liste bag
				addthissimu()
				duplicate/O simu_proba transindex
				makeindex/R simu_proba transindex
				// recup le delta PPM
				progressbarre[k][3]= labelpro[transindex[0]]
				// recup des indices pour chaques éléments
				progressbarre[k][4,118]= molecule[y-4]
			else
				molecule=0
				tuelesspec()
				genestringmol()
				addthissimu()
				progressbarre[k][3]=0
				progressbarre[k][4,118]= molecule[y-4]
			endif
			doupdate
			inco=gotonextintensity(inco)
			incox=wave0stripped(prochepic(inco,0))
			incoxdm=defmass(incox)
			incoy=wave1stripped(prochepic(inco,0))
	endfor
	inco=gotopreviousintensity(inco)
	incox=wave0stripped(prochepic(inco,0))
	incoxdm=defmass(incox)
	incoy=wave1stripped(prochepic(inco,0))
	
	
	/////////////////////////////////////////////////////////
	// Sauvegarde du molBag correspondant à l'attribution
	// et des chargeBag, listBag, cibleBag associés.
	//
	// On place le résultat dans le molbag current
	// Possibilité d'imaginer un systeme sauvegarde
	// différent.
	string nom="current"
	Majliste()
	Majlegende()
	
	/////////////////////////////////////////////////////////
	// Fabrique la wave de Rapport
	//
	// Fabrique la wave somme sur les collones pour detecter les éléments
	// présent ou non et fait le tri
	Make/O/N=118 somme
	Matrixop/O somme=sumcols(progressbarre)
	for(k=0, j=0;k<118;k+=1)
		if(somme[0][j]==0)
			deletepoints/M=1 j,1, progressbarre, somme, titres
		else
			j+=1
		endif
	endfor
	// Fabrique les titres
	Make/O/T/N=(1,118) titres
	Titres[0][0]="Proximite"
	Titres[0][1]="MassInconnue"
	Titres[0][2]="Intensitepic"
	Titres[0][3]="Delta"
	Titres[0][4,117]=elements[y-4]
	// Fait le Rapport
	Make/O/T/N=(n+1,j+1) Rapport
	Rapport[0][0]="Formule"
	Rapport[0][1,j+1]=titres[0][y-1]
	Rapport[1,n][1,j+1]=num2str(progressbarre[x-1][y-1])
	Rapport[1,n][0]=List_formules[numpnts(List_formules)-n+x-1]
	string masseprecise
	for(k=1;k<n+1;k+=1)
		sprintf masseprecise "%.6f", progressbarre[k-1][1]
		Rapport[k][2]=masseprecise
	endfor
	
	/////////////////////////////////////////////////////////
	// Demarage compare Attrib
	//
	// Kill la windows de barre de progression
	If(Wintype("showprogbar")!=0)
		KillWindow showprogbar
	endIf
	
	//  KILL la window de comparaison d'attribution,
	// Ceci entraine un ménage dans les waves.
	If(Wintype("AttributionComparaison")!=0)
		KillWindow AttributionComparaison
	endIf
	
	// Lance les fonctions d'initialisation de compareAttrib
	prepareCompareAttribution(nom)
	removeFakes()
	// Démarre l'interface GUI de compare Attrib
	interfaceCompareAttribution()

	killwaves somme, titres
end //fullauto(nprochains,seuil)

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compare Attrib
//
//	I.b Nouveaux éléments
//		Elements d'implémentation de la nouvelle interface de compare
//	Attrib. Rend obsolete les précédentes fonctions de visualisation
//	dévelopées par Roland (février 2014).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
//  prepareCompareAttribution(nomMolBag)
//	Fonction de préparation de l'environnement pour l'interface GUI compare
//	attrib en métant en place les waves qui vont bien.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function prepareCompareAttribution(nomMolBag)
	// Nom du mol bage que l'on veut avoir
	String nomMolBag
	// Variable globale conteneur du pointeur vers la donnée actuellment sélectionée
	Variable /G pointeurAttribCompareSelec = 0
	// Wave contenant les elements
	Wave /T elements
	////////////////////////////////////////////////////////////////////////////////
	//
	// CHARGEMENT DU MOLBAG DANS LE MODULE DE COMPARAISON
	//
	////////////////////////////////////////////////////////////////////////////////
	String nomWaveMolBag = "Molbag_" + nomMolBag
	String nomWaveCibleBag = "CibleBag_" + nomMolBag
	String nomWaveListBag = "ListBag_" + nomMolBag
	Wave /T waveListBag = $nomWaveListBag
	// Wave de couleurs 
	Make /O /N=14 rougeCompareAttrib={0,65000,0,0,33000,33000,60000,12500,0,0,46000,16000,33000,46000,65000}
	Make /O /N=20 vertCompareAttrib={0,0,33000,16000,0,33000,49000,53000,0,65000,48000,33000,33000,3000,16000}
	Make /O /N=20 bleuCompareAttrib={0,0,0,33000,33000,33000,5000,14000,65000,0,17000,33000,0,3000,0}
	// Recup du mol bag contenant les informations sur l'attribution
	Duplicate /O $nomWaveMolBag molbagCompareEnCours
	Duplicate /O $nomWaveCibleBag ciblebagCompareEnCours
	Duplicate /O waveListBag listbagCompareEnCours
	// Wave contenant les éléments dans le molbag
	Make /T /O /N=0 molbagCompareEnCoursElements
	Wave /T molbagCompareEnCoursElements
	Make /O /T /N=0 compareEnCoursElementsWaves
	nettoyageWavesElements()
	// Fabrique la wave somme sur les collones pour detecter les éléments présent ou non et fait le tri
	Matrixop/O somme=sumRows(molbagCompareEnCours)
	Variable k = 0
	Variable nombreElements = 0
	For(k=0 ; k<114 ; k+=1)
		If(somme[k]!=0)
			// Ajoutte point
			InsertPoints nombreElements,1,molbagCompareEnCoursElements,compareEnCoursElementsWaves
			nombreElements +=1
			// On a un element il faut construire les wave elements
			molbagCompareEnCoursElements[nombreElements-1] = elements[k]
			String nomWaveElement = "element_" + elements[k]
			Variable nombreMolecules = DimSize(molbagCompareEnCours,1)
			Make /O/D /N=(nombreMolecules) $nomWaveElement
			Wave waveElementEncours = $nomWaveElement
			// On récupère les données, repassage horizontal vers verticale, car le mol bag est la seule entité fonctionnant sur cette logique horizontale
			Variable j
			For(j=0 ; j< nombreMolecules ; j+=1)
				waveElementEncours[j] = molbagCompareEnCours[k][j]
			endFor
			compareEnCoursElementsWaves[nombreElements-1] = nomWaveElement
		Else
		endIf
	endfor
	Duplicate /O /R=[][0] ciblebagCompareEnCours Proximite
	Duplicate /O /R=[][1] ciblebagCompareEnCours MassInconnue
	Duplicate /O /R=[][2] ciblebagCompareEnCours Intensitepic,logintensitePic
	Duplicate /O /R=[][3] ciblebagCompareEnCours Delta
	logintensitePic=log(intensitePic)
End

////////////////////////////////////////////////////////////////////////////////////////////////////////
// miseAJourFenetreCompareAttrib()
//	Met à jour la fenêtre principale de compareAttrib
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function miseAJourFenetreCompareAttrib(nom)
	string nom
	Wave /T compareEnCoursElementsWaves
	Wave /T molbagCompareEnCoursElements
	Wave Delta
	Wave Proximite
	Wave MassInconnue
	Wave Intensitepic
	Wave logintensitePic
	
	Variable j
	For(j=0 ; j<numpnts(compareEnCoursElementsWaves) ; j+=1)
		RemoveFromGraph /W=AttributionComparaison#nombreAtome $compareEnCoursElementsWaves[j]
		String nomCheckBox = "affichageElement" + num2str(j)
		CheckBox $nomCheckBox disable=1
	endFor
	
	TitleBox selectionMolecule win=AttributionComparaison,title="Molécule : NONE"
	TitleBox affichageInfosMolecule1 win=AttributionComparaison,title="Masse : NONE\rIntensite : NONE\rProximite : NONE\rDelta : NONE"
	
	prepareCompareAttribution(nom)
	removeFakes()
	
	String textLengend = ""
	Wave rougeCompareAttrib
	Wave vertCompareAttrib
	Wave bleuCompareAttrib
	For(j=0 ; j<numpnts(compareEnCoursElementsWaves) ; j+=1)
		AppendToGraph /W=AttributionComparaison#nombreAtome $compareEnCoursElementsWaves[j] vs MassInconnue
		ModifyGraph /W=AttributionComparaison#nombreAtome mode[j]=3,zmrkSize[0]={logintensitePic,*,*,1,10},marker[j]=41,mrkThick[j]=1.5,rgb[j]=(rougeCompareAttrib[j],vertCompareAttrib[j],bleuCompareAttrib[j])
		String nomCheckBoxCreation = "affichageElement" + num2str(j)
		CheckBox $nomCheckBoxCreation  win=AttributionComparaison,title=molbagCompareEnCoursElements[j],size={35,15},pos={1055 + (j-mod(j,3))/3 *40 , 5+ (mod(j,3)) * 18},value=1,proc=affichageElementChanged,disable=0
		textLengend = textLengend + "\s(#" + num2str(j) + ") "+ molbagCompareEnCoursElements[j] + " \r"
	endFor
	ModifyGraph /W=AttributionComparaison#nombreAtome margin(top)=8, margin(left)=70, margin(right)=10,lblMargin(left)=30
	Label/W=AttributionComparaison#nombreAtome left "\\Z14Nb atomes"
	Label/W=AttributionComparaison#nombreAtome bottom "\\Z16Masse"
	
	textLengend = textLengend[0,strlen(textLengend)-3]
	Legend /C /W=AttributionComparaison#nombreAtome /N=legendNombreAtome textLengend
	
End //miseAJourFenetreCompareAttrib()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// hookAttributionCompare(s)
//	Fonction de gestion des évênements sur la fenêtre CompareAttrib.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function hookAttributionCompare(s)
	// strucutre contenant toutes les infos sur les évênements de la fenêtre
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	switch(s.eventCode)
		case 2 :
			// La fenetre va etre fermée va falloir faire du ménage
			// retrun0 pour permettre à Igor de gerer lui aussi l'evenement
			hookResult=0
			If(Wintype("GraphsAttributionComparaison")!=0)
				KillWindow GraphsAttributionComparaison
			endIf
			menageCompareAttrib()
			break
		case 3:
			// Click*
			// IMPORTANT return à 0 pour garantir la possibilité de double click pour modif graph
			hookResult=0
			traitementClickCompareAttrib(s.mouseLoc.h , s.mouseLoc.v)		
			If(s.eventMod==9)
				 supprimeMoleculeCompareAttrib("ctrlClick")
			endIf
			break
		case 11:
			// On a une entrée clavier
			If(s.keycode < 32 && s.keycode >27)
				// On gere les fleches nous meme, on ne laisse pas Igor le faire car il ne fait pas ce que l'on attend de lui
				If(s.keycode == 30 || s.keycode == 29)
					// On va de l'avant
					goToNextAttribCompare()
				Else
					// On va en arriere
					goToPreviousAttribCompare()
				endIf
				hookResult=1
			Else
				// Igor doit gerer l'evenement
				hookResult=0
			endIf
			break
	endswitch
	return hookResult
End //hookAttributionCompare(s)

////////////////////////////////////////////////////////////////////////////////////////////////////////
//  chargerMolbagAttribCompare(ctrlName,popNum,popStr) 
//	Chargement du molbag sélectionné dans le molbag_Current pour affichage
//	dans compareAttrib. Le précédent Current étant écrasé il faut sauvegarder
//	le précédent molbag avant si on veut le conserver.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function chargerMolbagAttribCompare(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	
	chargementMolBagAttribCompare(popStr)
End // chargerMolbagAttribCompare(ctrlName,popNum,popStr) 

////////////////////////////////////////////////////////////////////////////////////////////////////////
// saveMolBageCurrentAttribCompare(nomBouton)
//	Sauvegarde du molbag Current et des fonctions associées.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function saveMolBageCurrentAttribCompare(nomBouton) : ButtonControl
	String nomBouton
	
	NVar pointeurAttribCompareSelec
	pointeurAttribCompareSelec = 0
	
	String nomSauvegarde = "SaveMolbag"
	Prompt nomSauvegarde "Nom de sauvegarde du Molbag"
	DoPrompt "Nom du Molbag pour sauvegarde, s'il-vous-plaît.", nomSauvegarde
	
	If(V_flag==0)
		Wave Molbag_Current
		Wave Chargebag_Current
		Wave Listbag_Current
		Wave Ciblebag_Current
		
		String nomComplet = "Molbag_" + nomSauvegarde
		Duplicate /O Molbag_Current $nomComplet
	
		nomComplet = "Chargebag_" + nomSauvegarde
		Duplicate /O Chargebag_Current $nomComplet
		
		nomComplet = "Listbag_" + nomSauvegarde
		Duplicate /O Listbag_Current $nomComplet
		
		nomComplet = "Ciblebag_" + nomSauvegarde
		Duplicate /O Ciblebag_Current $nomComplet
	endIf
	
	String Camembert = recupListeMolBag()
	PopupMenu choixMolBag value=#camembert
	PopupMenu choixMolBag popmatch="Molbag_" + nomSauvegarde
	
	miseAJourFenetreCompareAttrib(nomSauvegarde)
End //saveMolBageCurrentAttribCompare(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// trierAttribCompare(nomMenu,indexSelect,strSelect)
//	Applique un tri sur les waves de l'actual molbag en fonction de la
//	wave sélectionée.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function trierAttribCompare(nomMenu,indexSelect,strSelect) : PopupMenuControl
	String nomMenu
	Variable indexSelect
	String strSelect
	
	NVar pointeurAttribCompareSelec
	
	Wave /T compareEnCoursElementsWaves
	Wave /T listbagCompareEnCours
	Wave Delta
	Wave Proximite
	Wave MassInconnue
	Wave Intensitepic 
	Wave logintensitePic
	
	Variable nombreElements = numpnts(compareEnCoursElementsWaves)
	Variable j = 0
	
	If(indexSelect < 5)
		// Pas Element
		Switch(indexSelect)
			case 1:
				// Tri en Masse
				For(j=0 ; j<nombreElements ; j+=1)
					Sort MassInconnue,$compareEnCoursElementsWaves[j]
				endFor
				Sort MassInconnue,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
			Break
			case 2:
				// Tri en Intensite
				For(j=0 ; j<nombreElements ; j+=1)
					Sort Intensitepic,$compareEnCoursElementsWaves[j]
				endFor
				Sort Intensitepic,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
			Break
			case 3:
				// Tri en deltaPPM
				For(j=0 ; j<nombreElements ; j+=1)
					Sort Delta,$compareEnCoursElementsWaves[j]
				endFor
				Sort Delta,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
			Break
			case 4:
				// Tri en Proximite
				For(j=0 ; j<nombreElements ; j+=1)
					Sort Proximite,$compareEnCoursElementsWaves[j]
				endFor
				Sort Proximite,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
			Break
		endSwitch
	Else
		// Element
		// Tri en Element 5-index
		Duplicate /O $compareEnCoursElementsWaves[indexSelect-5] elementTriOLD
		
		For(j=0 ; j<nombreElements ; j+=1)
			Sort elementTriOLD,$compareEnCoursElementsWaves[j]
		endFor
		Sort elementTriOLD,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
		
		KillWaves /Z elementTriOLD
	endIf
	
	pointeurAttribCompareSelec=0
	changeSelectionAttribCompare()
End // trierAttribCompare(nomMenu,indexSelect,strSelect)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// prevCompareAttribBut(nomBouton)
//	Apelle la fonction qui place les curseurs sur la molécule précédante en fonction
//	du critère contenu dans la wave selectionée pour le tri.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function prevCompareAttribBut(nomBouton) : ButtonControl
	String nomBouton
	goToPreviousAttribCompare()
End // prevCompareAttribBut(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// nextCompareAttribBut(nomBouton)
//	Apelle la fonction qui place les curseurs sur la molécule suivante en fonction
//	du critère contenu dans la wave selectionée pour le tri.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function nextCompareAttribBut(nomBouton) : ButtonControl
	String nomBouton
	goToNextAttribCompare()
End //nextCompareAttribBut(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// graphsCompareAttrib(nomBouton)
//	Appelle la fonction créant l'interdace d'affichage des graphs de
//	comparaison.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function graphsCompareAttrib(nomBouton) : ButtonControl
	String nomBouton
	interfaceGraphsCompAttrib()
End //graphsCompareAttrib(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// affichageElementChanged(nomCheckBox,checked)
//	Fonction appelée lors d'un changements des check box pilottant
//	l'affcihage des différents éléments.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function affichageElementChanged(nomCheckBox,checked) : CheckBoxControl
	String nomCheckBox
	Variable checked
	
	Variable extractionNumero = str2num(nomCheckBox[16,strlen(nomCheckBox)-1])
	
	If(checked == 0)
		ModifyGraph /W=AttributionComparaison#nombreAtome hideTrace[extractionNumero]=2
	Else
		ModifyGraph /W=AttributionComparaison#nombreAtome hideTrace[extractionNumero]=0
	EndIf
End //affichageElementChanged(nomCheckBox,checked)

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gestion Molbag
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
//  chargementMolBagAttribCompare(nomMolBagComplet)
//	Fonction de chargement du molbag nomMolBagComplet Molbag_blabla
//	et des fichiers associés dans les différentes waves de travail de compare attrib.
////////////////////////////////////////////////////////////////////////////////////////////////////////	
Function chargementMolBagAttribCompare(nomMolBag)
	String nomMolBag
		
	PopupMenu choixTri popmatch="Proximité"
	
	rechargerCibleBag(nomMolBag)
	
	miseAJourFenetreCompareAttrib(nomMolBag)
	
	If(Wintype("GraphsAttributionComparaison")!=0)
		KillWindow GraphsAttributionComparaison
	endIf
End //chargementMolBagAttribCompare(nomMolBagComplet)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// sauvegardeMolBagApresAttrib(nom)
//	Sauvegarde du molbag Current sous le nom "nom".
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function sauvegardeMolBagApresAttrib(nom)
	String nom
	
	// Wave contenant toutes les données de la dernière attribution
	Wave resultatAttribution = progressBarre
	
	String nomDuMolBag = "Molbag_" + nom
	Variable nombreMolecules = DimSize(resultatAttribution,0)
	Make /O/D/N=(114,nombreMolecules) $nomDuMolBag
	Wave molBagEnConstruction = $nomDuMolBag
	
	// On récupère les données
	Variable k
	For(k=0 ; k < nombreMolecules ; k+=1)
		Variable j
		For(j=0 ; j< 114 ; j+=1)
			molBagEnConstruction[j][k] = resultatAttribution[k][j+4]
		endFor
	endFor
End //sauvegardeMolBagApresAttrib(nom)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// sauvegardeChargeBagApresAttrib(nom)
//	Sauvegarde du chargebag Current sous le nom "nom" avec la fonction
//	de Maylis.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function sauvegardeChargeBagApresAttrib(nom)
	String nom
	// Recupere la charge avec fonction Maylis
	chopelacharge(nom)
End //sauvegardeChargeBagApresAttrib

////////////////////////////////////////////////////////////////////////////////////////////////////////
// sauvegardeListBagApresAttrib(nom)
//	Sauvegarde du listbag Current sous le nom "nom" avec la fonction
//	de Maylis.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function sauvegardeListBagApresAttrib(nom)
	String nom
	wave/T list_formules

	String nomdelaliste
	nomdelaliste="Listbag_"+nom
	Duplicate/O list_formules $nomdelaliste
End //sauvegardeListBagApresAttrib

////////////////////////////////////////////////////////////////////////////////////////////////////////
// sauvegardeCibleBagApresAttrib(nom)
//	Sauvegarde du ciblebag Current sous le nom "nom" en piochant
//	dans la wave progressBarre issue de l'attribution.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function sauvegardeCibleBagApresAttrib(nom)
	String nom
	
	// Wave contenant toutes les données de la dernière attribution
	Wave resultatAttribution = progressBarre
	
	String nomDuCibleBag = "Ciblebag_" + nom
	Variable nombreMolecules = DimSize(resultatAttribution,0)
	Make /O/D/N=(4,nombreMolecules) $nomDuCibleBag
	Wave cibleBagEnConstruction = $nomDuCibleBag
	
	Duplicate /R=[][0,4] /O resultatAttribution, cibleBagEnConstruction
End //sauvegardeCibleBagApresAttrib

////////////////////////////////////////////////////////////////////////////////////////////////////////
// rechargerCibleBag()
//	Charge le cibleBAg contenant les infos masse, intensite, proximite et PPM
//	dans les waves correspondante pour utilisation dans compare attrib.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function rechargerCibleBag(nom)
	string nom
	string sciblebag="Ciblebag_"+nom
	Wave Proximite
	Wave MassInconnue
	Wave Intensitepic
	Wave Delta
	
	Duplicate /O /R=[][0] $sciblebag Proximite
	Duplicate /O /R=[][1] $sciblebag MassInconnue
	Duplicate /O /R=[][2] $sciblebag Intensitepic
	Duplicate /O /R=[][3] $sciblebag Delta
End // rechargerCibleBag()


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Déplacement dans données
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// goToNextAttribCompare()
//	Place les curseurs sur la molécule suivante en fonction
//	du critère contenu dans la wave selectionée pour le tri.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function goToNextAttribCompare()
	NVar pointeurAttribCompareSelec
	
	Wave Massinconnue
	If(pointeurAttribCompareSelec<numpnts(Massinconnue)-1)
		pointeurAttribCompareSelec+=1
		changeSelectionAttribCompare()
	endIf
End // goToNextAttribCompare()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// goToPreviousAttribCompare()
//	Place les curseurs sur la molécule précédante en fonction
//	du critère contenu dans la wave selectionée pour le tri.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function goToPreviousAttribCompare()
	NVar pointeurAttribCompareSelec
	
	If(pointeurAttribCompareSelec>0)
		pointeurAttribCompareSelec-=1
		changeSelectionAttribCompare()
	endIf
End // goToreviousAttribCompare()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// changeSelectionAttribCompare()
//	Place les curseurs sur la position contenue dans la variable globale
//	pointeurAttribCompareSelec.
//	Actualise les curseurs de l'interface compare attrib mais aussi
//	du curseur dans le spectre, dans le dmvm et enfin dans l'interface
//	graphsCompareAttrib si celle-ci existe.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function changeSelectionAttribCompare()
	NVar pointeurAttribCompareSelec
	
	SVAR nomWaveListBag = nomWaveListBag
	Wave /T listbagCompareEnCours = listbagCompareEnCours
	
	Cursor /W=AttributionComparaison#Proximite /S=0 /C=(45000,0,0) curs1 Proximite pointeurAttribCompareSelec
	Cursor /W=AttributionComparaison#DeltaMass /S=0 /C=(45000,0,0) curs2 Delta pointeurAttribCompareSelec
	
	// Les noms des curseurs sont reglementés dans Igor pour en avoir plusieurs sur le même graph. Me demandez pas pourquoi, surement une norme européenne...
	// En tout cas, cela permet de réviser son alphabet - toujours utile surtout en cas de voyage linguiste en Patagonie Orientale.
	// En effet, dans ce pays l'ignorance de son alphabet peut être condamné par la peine de mort. L'execution de la sentece prenant la forme d'un combat ritualisé
	// face à une horde de pécaris enragés.
	string reserveMemoireCursorsName = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	
	// Placement des curseurs sur les waves elements
	// wave des noms des éléments présents
	Wave /T compareEnCoursElementsWaves
	Variable o
	For(o=0 ; o<numpnts(compareEnCoursElementsWaves) ; o+=1)
		String nomCursor = reserveMemoireCursorsName[o]
		Cursor /W=AttributionComparaison#nombreAtome /S=2 /C=(45000,0,0) $nomCursor $compareEnCoursElementsWaves[o] pointeurAttribCompareSelec
	endFor
	
	Wave Delta
	Wave Proximite
	Wave MassInconnue
	Wave Intensitepic
	// Change nom de la molecule
	TitleBox selectionMolecule win=AttributionComparaison,title=("Molécule : " + listbagCompareEnCours[pointeurAttribCompareSelec])
	TitleBox affichageInfosMolecule1 win=AttributionComparaison,title="Masse : "+ num2str(MassInconnue[pointeurAttribCompareSelec]) + "\rIntensite : " + num2str(Intensitepic[pointeurAttribCompareSelec]) + "\rProximite : " + num2str(Proximite[pointeurAttribCompareSelec]) + "\rDelta : " + num2str(Delta[pointeurAttribCompareSelec])
	
	// On de place les Curseurs agregateur et dmvm
	Wave incox
	Wave incoy
	Wave incoxdm
	incox=MassInconnue[pointeurAttribCompareSelec]
	incoxdm=defmass(incox)
	incoy=Intensitepic[pointeurAttribCompareSelec]
	
	// Si la fenêtre de Graphs compare attrib existe alors on déplace ses curseurs
	If(winType("GraphsAttributionComparaison") != 0)
		placerCursorsCompareGraphs()
	endIf
End // changeSelectionAttribCompare()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// traitementClickCompareAttrib(positionX,positionY)
//	Detecte la masse correspondant au point le plus proche de la zone
//	de click passée en paramètres.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function traitementClickCompareAttrib(positionX,positionY)
	Variable positionX
	Variable positionY
	
	NVar pointeurAttribCompareSelec
	
	pointeurAttribCompareSelec = str2num(stringbykey("HITPOINT",TraceFromPixel (positionX,positionY,"")))

	changeSelectionAttribCompare()
End //traitementClickCompareAttrib(positionX,positionY)


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// utilitaires CompareAttrib
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// genereListeFromTextWave(textWave)
//	Genere la liste pour le popup du menu de tri des données à partir
//	d'une wave texte.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function /S genereListeFromTextWave(textWave)
	Wave /T textWave
	
	String resultat = "\""
	resultat = resultat + "Masse;Intensité;PPM;Proximité;"
	Variable j
	For(j=0 ; j<numpnts(textWave) ; j+=1)
		resultat = resultat + textWave[j] + ";"
	endFor
	
	resultat = resultat + "\""
	return resultat
End // genereListeFromTextWave(textWave)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// recupListeMolBag()
//	Grâce a un ours au yeux qui tirent des laser le Canada nous récupère
//	les différents molbag présent en mémoire
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function /S recupListeMolBag()
	String OCanada = "\"" + WaveList("Molbag_*",";","") + "\""
	return OCanada
End // recupListeMolBag()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// nettoyageWavesElements()
//	Nettoyage réccursif des waves contenant les éléments.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function nettoyageWavesElements()
	String elements = wavelist("element_*", ";", "")
	variable n=itemsinlist(elements), k
	string lawave
	for(k=0;k<n;k+=1)
		lawave=stringfromlist(k, elements)
		killwaves/Z $lawave
	endfor
End // nettoyageWavesElements()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// menageCompareAttrib()
//	Nettoyage des waves de compare attrib.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function menageCompareAttrib()
	Killwindow AttributionComparaison#nombreAtome
	Killwindow AttributionComparaison#DeltaMass
	Killwindow AttributionComparaison#Proximite
	Killwaves /Z ciblebagCompareEnCours,molbagCompareEnCoursElements,compareEnCoursElementsWaves,molbagCompareEnCours,listbagCompareEnCours
	Killwaves /Z Delta,Proximite,MassInconnue,Intensitepic,logintensitepic,rougeCompareAttrib, vertCompareAttrib, bleuCompareAttrib
	nettoyageWavesElements()
End // menageCompareAttrib()


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interface affichage Graphs
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// hookAttributionCompareGraphs(s)
//	Gestion des événements de la fenetre d'interface d'affichage des Graphs
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function hookAttributionCompareGraphs(s)
STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	
	// Récup la hauteur de la window et sa largueur
	Variable hauteur = abs(s.winRect.top - s.winRect.bottom)
	Variable largueur = abs(s.winRect.left - s.winRect.right)
	
	Variable hauteurDispo = hauteur -40
	Variable hauteurGraph = hauteurDispo / 2
	
	Variable largueurDispo = largueur-50
	Variable largueurGraph = largueurDispo / 4
	
	switch(s.eventCode)
		case 2 :
			// La fenetre va etre fermée il faut faire du nettoyage, et réactiver l'option de tri
			// dans l'interface principale de compare Attrib.
			// retrun0 pour permettre à Igor de gerer aussi l'evenement
			hookResult=0
			PopupMenu choixTri win=AttributionComparaison,disable=0
			 nettoyageCompareAttribGraphs()
			break
		case 3:
			// Click*
			// IMPORTANT return à 0 pour garantir la possibilité de double click pour modif graph
			hookResult=0
			// Par les tests successifs sur la position de la souris on detecte si le click est sur un graph à selection (utilisant un jeu partiel de données)
			// ou sur un graph à jeu de données complet, important pour identification de la molécule sélectionée.
			If(s.mouseLoc.v > 35 && s.mouseLoc.h>45)	
				hauteurGraph = abs(s.winRect.top - s.winRect.bottom)
				largueurGraph = abs(s.winRect.left - s.winRect.right)
					
				If(s.winRect.left >( 25 + 3*largueurGraph))
					traitementClickCompareGraphs(s.mouseLoc.h , s.mouseLoc.v,1)
				Else
					If(s.winRect.left > 25 + 2*largueurGraph && s.winRect.top < 40)
						traitementClickCompareGraphs(s.mouseLoc.h , s.mouseLoc.v,1)
					Else
						traitementClickCompareGraphs(s.mouseLoc.h , s.mouseLoc.v,0)
					endIf
				endIf
			endIf	
			break
		case 6:
			// Redimensionne la fenetre
			// on recalcule les dimensions des différents graphiques
			MoveSubwindow /W=GraphsAttributionComparaison#graph11 fnum=(45+0*largueurGraph , 35+0*hauteurGraph , 45+1*largueurGraph , 35+1*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph21 fnum=(45+0*largueurGraph , 35+1*hauteurGraph , 45+1*largueurGraph , 35+2*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph12 fnum=(45+1*largueurGraph , 35+0*hauteurGraph , 45+2*largueurGraph , 35+1*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph22 fnum=(45+1*largueurGraph , 35+1*hauteurGraph , 45+2*largueurGraph , 35+2*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph13 fnum=(45+2*largueurGraph , 35+0*hauteurGraph , 45+3*largueurGraph , 35+1*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph23 fnum=(45+2*largueurGraph , 35+1*hauteurGraph , 45+3*largueurGraph , 35+2*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph14 fnum=(45+3*largueurGraph , 35+0*hauteurGraph , 45+4*largueurGraph , 35+1*hauteurGraph)
			MoveSubwindow /W=GraphsAttributionComparaison#graph24 fnum=(45+3*largueurGraph , 35+1*hauteurGraph , 45+4*largueurGraph , 35+2*hauteurGraph)
			
			break
	endswitch
	return hookResult	
End // hookAttributionCompareGraphs(s)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// selectionNbOGraphsCompareAttrib(nomSlider, valeur, event)
//	Slider permettant la selection dans des waves de selection des
//	molécules remplissant la condition sur le nombre d'oxygène.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function selectionNbOGraphsCompareAttrib(nomSlider, valeur, event) : SliderControl
		String nomSlider
		Variable valeur
		Variable event
		
		// On rafrachit le title box servant à masque le -1 pour infos selection tous Oxygène
		// Sinon par le fait qu'il y a eu un evenement le slider repasse par dessus...
		TitleBox legendSlider win=GraphsAttributionComparaison
		
		selectionFonctionNbO()
End // selectionNbOGraphsCompareAttrib(nomSlider, valeur, event)

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Recup' dimension de l'écran
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
////////////////////////////////////////////////////////////////////////////////////////////////////////
// extractionHauteurEcran()
//	Recup de la hauteur de l'écran principal
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function /D extractionHauteurEcran()
	
	String InfosExtraite = StringByKey("SCREEN1",IgorInfo(0))
	String rectExtract = InfosExtraite[strsearch(InfosExtraite,"=",10)+1,strlen(InfosExtraite)-1]
	
	Variable largueur = str2num(rectExtract[findPosStringVirgule(rectExtract,3)+1 , strlen(rectExtract)-1]) - str2num(rectExtract[findPosStringVirgule(rectExtract,1)+1 , findPosStringVirgule(rectExtract,2)-1])
	return largueur
End //extractionHauteurEcran()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// extractionLargueurEcran()
//	Recup de la largueur de l'écran principal
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function /D extractionLargueurEcran()
	
	String InfosExtraite = StringByKey("SCREEN1",IgorInfo(0))
	String rectExtract = InfosExtraite[strsearch(InfosExtraite,"=",10)+1,strlen(InfosExtraite)-1]
	
	Variable hauteur = str2num(rectExtract[findPosStringVirgule(rectExtract,2)+1 , findPosStringVirgule(rectExtract,3)-1]) - str2num(rectExtract[0 , findPosStringVirgule(rectExtract,1)-1])
	return hauteur
End //extractionLargueurEcran()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// findPosStringVirgule(stringRecherche,virguleVoulue)
//	Extraction infos dans un string composé d'une liste séparée par
//	des virgules.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function findPosStringVirgule(stringRecherche,virguleVoulue)
	String stringRecherche
	Variable virguleVoulue
	
	Variable nbVirgulesTrouvees = 0
	
	Variable memoirePositionRecherche = 0
	Do
		memoirePositionRecherche=strsearch(stringRecherche,",",memoirePositionRecherche+1)
		nbVirgulesTrouvees +=1
	While(nbVirgulesTrouvees<virguleVoulue)
	
	Return memoirePositionRecherche
End //findPosStringVirgule(stringRecherche,virguleVoulue)


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Deplacement dans donnes Graphs
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// traitementClickCompareGraphs(positionX,positionY,mode)
//	Gestion du click dans la zone des graphs. Detecte la position de
//	la donnée la plus proche et en déduit une position dans les donnees
//	en rensignant la variable globale pointeurAttribCompareSelec. Cette
//	variable contient l'index de la molécule sélectionée dans le jeu de donnée
//	complet.
//	
//	La variable mode permet à cette fonction de savoir si le click a eu lieu
//	sur une graphs représentant une selection ou si le graphs ou a lieu le
//	click correspond à un jeu de données complète.
//	Permet le calcul des positions de façon cohérente pour propagation dans
//	les différents graphs et interfaces (vers comapre attrib).
///////////////////////////////////////////////////////////////////////////////////////////////////////
Function traitementClickCompareGraphs(positionX,positionY,mode)
	Variable positionX
	Variable positionY
	Variable mode //=0 click dans zone selection_ ; 1 zone sans selection en O
	
	NVar pointeurAttribCompareSelec
	
	// Pointeur dans la selection
	Variable pointeurGraphs
	pointeurGraphs = str2num(stringbykey("HITPOINT",TraceFromPixel (positionX,positionY,"")))
	
	Wave selection_masse
	Wave MassInconnue
	
	Variable masseRecherchee
	
	If(mode==0)		
		// on transforme de la selection vers le jeu de données entier
		masseRecherchee = selection_masse[pointeurGraphs]
		
		FindValue /V=(masseRecherchee) MassInconnue
		pointeurAttribCompareSelec = V_value
	Else		
		// Pointeur dans compareAttrib
		pointeurAttribCompareSelec = pointeurGraphs
	endIf
		
	// On récupère le nom de la molécule
	Wave /T listbagCompareEnCours
	TitleBox selectionMoleculeGraphs win=GraphsAttributionComparaison, title=("Molécule : " + listbagCompareEnCours[pointeurAttribCompareSelec])
	
	// On bouge le cursor dans compare atrib et agregateur
	changeSelectionAttribCompare()	
End //traitementClickCompareGraphs(positionX,positionY,mode)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// placerCursorsCompareGraphs()
//	Place les curseurs de l'interface Graphs en fonction de la variable
//	globale pointeurAttribCompareSelec. Fait une conversion pour les graphs
//	avec plus de points.
//////////////////////////////////////////////////////////////////////////////////////////////////////// 
Function placerCursorsCompareGraphs()
	Wave MassInconnue
	Wave selection_masse

	// position dans les données générales
	NVar pointeurAttribCompareSelec
	// on va recherhcée la masse correspondante dans le jeu de données de selection.
	Variable masseRecherchee = MassInconnue[pointeurAttribCompareSelec]
	FindValue /V=(masseRecherchee) selection_masse
	Variable pointeurGraphsWaveSelection_ = V_value
	
	If(V_value==-1)
		// La molécule sélectionée n'est pas dans la selection -> on ne doit pas afficher de curseurs sur les graphs de selection
		Cursor /K /W=GraphsAttributionComparaison#graph11 cursGraphs11
		Cursor /K /W=GraphsAttributionComparaison#graph21 cursGraphs21
		Cursor /K /W=GraphsAttributionComparaison#graph12 cursGraphs12
		Cursor /K /W=GraphsAttributionComparaison#graph22 cursGraphs22
		Cursor /K /W=GraphsAttributionComparaison#graph23 cursGraphs23
	Else
		// la molécule est dans les sélection, placement des curseurs dans les graphs de selection.
		Cursor /W=GraphsAttributionComparaison#graph11 /S=2 /C=(60000,60000,0) cursGraphs11 selection_HsurC pointeurGraphsWaveSelection_
		Cursor /W=GraphsAttributionComparaison#graph21 /S=2 /C=(60000,60000,0) cursGraphs21 selection_NsurC pointeurGraphsWaveSelection_
		Cursor /W=GraphsAttributionComparaison#graph12 /S=2 /C=(60000,60000,0) cursGraphs12 selection_DBE pointeurGraphsWaveSelection_
		Cursor /W=GraphsAttributionComparaison#graph22 /S=2 /C=(60000,60000,0) cursGraphs22 selection_composite pointeurGraphsWaveSelection_
		Cursor /W=GraphsAttributionComparaison#graph23 /S=2 /C=(60000,60000,0) cursGraphs23 selection_excesMethylene pointeurGraphsWaveSelection_
	endIf	
	
	// graphs sur jeu de données entier, affiche dans tous les cas
	Cursor /W=GraphsAttributionComparaison#graph13 /S=2 /C=(60000,60000,0) cursGraphs13 Zubarev pointeurAttribCompareSelec
	Cursor /W=GraphsAttributionComparaison#graph14 /S=2 /C=(60000,60000,0) cursGraphs14 VanKrev_HsurC pointeurAttribCompareSelec
	Cursor /W=GraphsAttributionComparaison#graph24 /S=2 /C=(60000,60000,0) cursGraphs24 VanKrev_HsurC pointeurAttribCompareSelec
End //placerCursorsCompareGraphs()


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Autres
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// nettoyageCompareAttribGraphs()
//	Nettoyage des waves de Graphs
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function nettoyageCompareAttribGraphs()
	KillWindow GraphsAttributionComparaison#graph11
	KillWindow GraphsAttributionComparaison#graph21
	KillWindow GraphsAttributionComparaison#graph12
	KillWindow GraphsAttributionComparaison#graph22
	KillWindow GraphsAttributionComparaison#graph13
	KillWindow GraphsAttributionComparaison#graph23
	KillWindow GraphsAttributionComparaison#graph14
	KillWindow GraphsAttributionComparaison#graph24
	
	// on vire les waves de selection en oxygène
	String selections = wavelist("selection_*", ",", "")
	Execute "killwaves " + selections
	
	Killwaves Zubarev, VanKrev_OsurC, VanKrev_HsurC, VanKrev_NsurC,VanKrevMass
End //nettoyageCompareAttribGraphs()




////revue le 30 octobre 2014 pour ajouter ADD CURRENT TARGET TO ROI
function getagreg(s)
	STRUCT WMWinHookStruct &s
variable hookres=0
nvar findindata,inco, treedepth, charge
wave/Z incox, incoxdm, cursint, incoy, wave1stripped
variable/g lastx, lasty, zoomxfactor=0, testcursorin=0
wave PeakSet_Current
svar hightlightedPeak
string locname="Peak_"+hightlightedPeak+"_Y"
make/o/d/n=1 pozmasse, pozint
pozmasse=axisvalfrompixel("agregateur","bottom",s.mouseLoc.h)
pozint=axisvalfrompixel("agregateur","left",s.mouseLoc.v)
if(s.eventCode==8)
return 0
endif
if(s.eventCode==3)
hookres=1
	if(s.eventMod==16 || s.eventMod==17)
	wave PeakSet_Current
		if(stringmatch(hightlightedPeak,"")||numpnts(PeakSet_Current)==0)
			if(itemsinlist(winlist("PeakQualCont",";","")))
				PopupContextualMenu "Calculate;Clear ALL Peaks;Hide Peak Window;Add Formula Peaks;Add a Peak for each ROI point;Add current m/z to ROI"
			else
				PopupContextualMenu "Calculate;Clear ALL Peaks;Show Peak Window;Add Formula Peaks;Add a Peak for each ROI point;Add current m/z to ROI"
			endif
		else
			PopupContextualMenu "Calculate;Clear ALL Peaks;Clear this peak;Fit this peak;Calculate Formula"
		endif
		strswitch(S_selection)	
			case "Calculate":
				Calculate(inco,incoy[0],charge)
				break
			case "Clear ALL Peaks":
				killALLpeaks(PeakSet_Current)
				break
			case "Clear this peak":
				KillThisPeak(PeakSet_Current,hightlightOnePeak(pozint[0], pozmasse[0], PeakSet_Current))
				break
			case "Fit this peak":
				FitThisPeak(PeakSet_Current,hightlightOnePeak(pozint[0], pozmasse[0], PeakSet_Current),hold2cons())
				break
			case "Show Peak Window":
				Execute "PeakQualCont()"
				break
			case "Hide Peak Window":
				Killwindow PeakQualCont
				break
			case "Calculate Formula":
				calculatefrompeak()
				break
			case "Add Formula Peaks":
				AddPeakFromCurrentSimu()
				break
			case "Add a Peak for each ROI point":
				AddPeakFromROI()
				break
			case "Add current m/z to ROI":
				AddCurrent2ROI()
				break
		endswitch
		return 1
	endif
	return 0
return hookres
endif
if(s.eventCode==4 && numpnts(cursint)==1)
//majcursint(pozint)
//dowindow/F agregateur
if (numpnts(PeakSet_Current)>0)
inco=hightlightOnePeak(pozint[0], pozmasse[0], PeakSet_Current)
incox=inco
incoxdm=defmass(incox)
	if(findindata>0)
		incoy=wave1stripped[mass2pnt(incox)]
	elseif(findindata==0)
		incoy=wavemax(wave1stripped)/2
	endif
endif
killcontrol sliderzoom
lastx=s.mouseLoc.h
lasty=s.mouseLoc.v
return 0
endif
if(s.eventCode==5)
getaxis/W=agregateur/Q left
testcursorin = pozint[0]<V_max && pozint[0]>V_min
getaxis/W=agregateur/Q bottom
testcursorin *= pozmasse[0]<V_max && pozmasse[0]>V_min
if(testcursorin)
	if(findindata>0)
		findclosepntagreg(pozmasse[0],pozint[0])
		majcursint(incoy[0])
		controlinfo/W=dmvm popup0
		generatetree(S_value,treedepth,inco)
	elseif(findindata==0)
		inco=pozmasse[0]
		incox=inco
		incoxdm=defmass(incox)
	endif
endif
	calcinfo()
endif
if(s.eventCode==22)
if (numpnts(wave1stripped)<200000)
	zoomHook(s.winName,"bottom",pozmasse[0],0.1*s.wheelDy,wavemin(wave0stripped),wavemax(wave0stripped))
	setaxis/W=$s.winName/A=2 left
	if(wintype("PeakQualCont"))
		getaxis/Q/W=$s.winName bottom
		setaxis/W=PeakQualCont bottom V_min, V_max
	endif
	if(s.eventMod==8)
		getaxis/Q/W=$s.winName bottom
		setaxis/W=dmvm bottom V_min, V_max
		setaxis/W=dmvm/A=2 left
	endif
	topbright()
	calcinfo()
else
print "Too many points, scrolling unavailable"
endif
return 0
endif
return 0
end

function AddCurrent2ROI()
wave roi0, roi1, roi2, incox, incoy, incoxdm
insertpoints 0,1,roi0, roi1, roi2
roi0[0]=incox[0]
roi1[0]=incoy[0]
roi2[0]=incoxdm[0]
sort roi0, roi0, roi1, roi2
calcinfo()
end

function cutROIfromagreg()
wave roi0,roi1,roi2
getmarquee/W=agregateur bottom
findlevel/Q roi0, max(V_left,wavemin(roi0))
variable b1=ceil(V_levelX)
findlevel/Q roi0, min(V_right, wavemax(roi0))
variable b2=floor(V_levelX)
variable n=b2-b1+1
if(stringmatch(S_marqueeWin,"agregateur")&& n>-1)
	deletepoints b1,n, roi0, roi1, roi2
	//print V_left,b1, V_right, b2, n
else
print n
endif
calcinfo()
end

function fullsavebag(nom)
string nom
variable t=ticks
singlemol2matmolbag(nom)
chopelacharge(nom)
sauvelaliste(nom)
sauveciblebag(nom)
//print "temps de sauvegarde du molbag de l'agregateur :", ticks-t
end

function sauveciblebag(nom)
string nom
wave List_targets
string nomdelaliste
nomdelaliste="Ciblebag_"+nom
duplicate/O List_targets $nomdelaliste
end

Function Ajouter(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			variable/g affichageAttribONOFF = 1
			addthissimu()
			majliste()
			majlegende()
			break
	endswitch

	return 0
End

function ManyPntsCal()
wave list_targets, index_formules, error1, error0, color
wave/T list_formules
variable n=numpnts(index_formules)
nvar degreedelm
duplicate/O index_formules, pntserror, pntsmes
pntserror=list_targets[x][1]-index_formules
pntsmes=list_targets[x][1]
colorit()
sort pntsmes, pntsmes, pntserror, color
make/O/T/N=(n,5) calcontrol=""
calcontrol[][0]=list_formules[x]
calcontrol[][1]=num2str(index_formules[x])
calcontrol[][2]=num2str(list_targets[x][1])
calcontrol[][3]=num2str(list_targets[x][3])
CurveFit/Q/M=2/W=0 poly (degreedelm), pntserror/X=pntsmes/D
error1=poly(W_coef,error0[x])
end

Function CalFromMolbag(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
wave error1, error0, color
nvar degreedelm
Duplicate/O molbag2fullindx(bag), pntserror, pntsmes
pntserror=ciblebag[x][1]-target2theo(ciblebag[x][1],ciblebag[x][3])
pntsmes=ciblebag[x][1]
sort pntsmes, pntsmes, pntserror
make/O/T/N=(numpnts(chargebag),5) calcontrol=""
calcontrol[][0]=listbag[x]
calcontrol[][1]=num2str(target2theo(ciblebag[x][1],ciblebag[x][3]))
calcontrol[][2]=num2str(ciblebag[x][1])
calcontrol[][3]=num2str(ciblebag[x][3])
CurveFit/Q/M=2/W=0 poly (degreedelm), pntserror/X=pntsmes/D
error1=poly(W_coef,error0[x])
end

Function SetVarProc_3(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			nvar massgap, calwidth
			delM(massgap,calwidth)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_7(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			ManyPntsCal()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_4(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave pointconf
nvar intenscumthres
	switch( ba.eventCode )
		case 2: // mouse up
			correctmass()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Window datamanager() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1507,595,1704,776)
	ListBox list0,pos={1,1},size={192,128}
	ListBox list0,help={"Displays available datasets. Click to change and refresh top ten peaks."}
	ListBox list0,listWave=root:liste_donnees,mode= 1,selRow= -1
	Button buttonLoadData,pos={3,130},size={59,25},proc=ButtonLoader,title="Load"
	Button buttonLoadData,help={"Load an XY text file with mass and intensities."}
	Button buttonLoadData,fColor=(32768,65280,0)
	Button buttonKiller,pos={102,155},size={35,25},proc=ButtonKill,title="Kill"
	Button buttonKiller,help={"Removes the selected sample data."}
	Button buttonKiller,fColor=(65280,0,0)
	Button buttonRename,pos={102,130},size={35,25},proc=ButtonDupeRname,title="Dupe"
	Button buttonRename,help={"Generates an other sample data from the current data."}
	Button buttonRename,fColor=(32768,65280,0)
	Button buttonFatNoise1,pos={4,157},size={59,25},proc=ActiveLaDataButton,title="Activate"
	Button buttonFatNoise1,help={"Remove the less intense points so the new intensity distribution is strictly convex"}
	Button buttonFatNoise1,fColor=(65280,43520,0)
	Button buttonLoadinROI,pos={64,135},size={34,38},proc=loadInRoi,title="To ROI"
	Button buttonLoadinROI,help={"Load an XY text file with mass and intensities."}
	Button buttonLoadinROI,fColor=(32768,65280,0)
EndMacro

Window isotopefinder() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(720,138.5,925.5,315.5) histoiso
	AppendToGraph/T histoiso
	AppendToGraph/B=top2 histoiso
	AppendToGraph/B=top3 histoiso
	ModifyGraph userticks(top)={defmassreliso12,elements}
	ModifyGraph userticks(top2)={defmassreliso13,elements}
	ModifyGraph userticks(top3)={defmassreliso23,elements}
	ModifyGraph mode=5
	ModifyGraph rgb=(0,0,0)
	ModifyGraph grid(top)=1,grid(top2)=1,grid(top3)=1
	ModifyGraph tick(left)=3
	ModifyGraph noLabel(left)=2
	ModifyGraph fStyle(top)=1,fStyle(top2)=1,fStyle(top3)=1
	ModifyGraph lblMargin(top)=32
	ModifyGraph axOffset(left)=-6,axOffset(bottom)=-1,axOffset(top)=0.5
	ModifyGraph axThick(left)=0
	ModifyGraph gridRGB(top2)=(0,52224,0),gridRGB(top3)=(65280,0,0)
	ModifyGraph axRGB(top)=(24576,24576,65280),axRGB(top2)=(0,52224,0),axRGB(top3)=(65280,0,0)
	ModifyGraph tlblRGB(top)=(24576,24576,65280),tlblRGB(top2)=(0,52224,0),tlblRGB(top3)=(65280,0,0)
	ModifyGraph alblRGB(top)=(24576,24576,65280),alblRGB(top2)=(0,52224,0),alblRGB(top3)=(65280,0,0)
	ModifyGraph lblPos(bottom)=25
	ModifyGraph freePos(top2)=-109
	ModifyGraph freePos(top3)=-81
	Label top "\\Z061\\Sst\\M\\Z06 - 2\\Snd"
	Label top2 "\\Z061\\Sst\\M\\Z06 - 3\\Srd"
	Label top3 "\\Z062\\Snd\\M\\Z06 - 3\\Srd"
	Button button0,pos={185,3},size={50,20},proc=ButtonZ2,title="Zoom2"
	Button button1,pos={57,3},size={50,20},proc=ButtonZ1,title="Zoom1"
	Button button2,pos={121,3},size={50,20},proc=ButtonZ0,title="Reset"
EndMacro

Window molmanager() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(267,57,421,685)
	ListBox list0,pos={1,31},size={152,371},proc=ListBoxProcAgregToElab
	ListBox list0,help={"Displays the list of formulas shown in the agregator. Click to select one."}
	ListBox list0,labelBack=(43520,43520,43520),fSize=12,frame=4,fStyle=0
	ListBox list0,listWave=root:listboxformules,selWave=root:selboxformules
	ListBox list0,titleWave=root:listboxformulestitles,mode= 1,selRow= -1
	ListBox list0,special= {0,18,0},widths={20,69,29,44},userColumnResize= 1
	Button Delete,pos={2,486},size={76,25},proc=DeleteMol,title="Erase this"
	Button Delete,help={"Removes selected formula from list and from agregator window"}
	Button Delete,fColor=(65280,0,0)
	Button centre,pos={2,436},size={76,25},proc=centrevue,title="Focus"
	Button centre,help={"Changes the agregator window settings to show the entire isotopical simulation. Focus when none selected restore default settings."}
	Button centre,fColor=(32768,65280,0)
	Button centre1,pos={77,436},size={76,25},proc=Switchmass,title="Iso Pattern"
	Button centre1,help={"Spans each peak of an isotopical simulation in the agregator window."}
	Button centre1,fColor=(32768,65280,0)
	Button autoscale,pos={2,461},size={76,25},proc=rescaleselected,title="Rescale this"
	Button autoscale,help={"Sets the maximum of selected isotopical simulation in the agregator equal to the height of the closest peak in the dataset"}
	Button autoscale,fColor=(65280,43520,0)
	Button rescaleALL,pos={77,461},size={76,25},proc=rescaleall,title="Rescale all"
	Button rescaleALL,help={"Sets the maximum of each isotopical simulation in the agregator equal to the height of the closest peak in the dataset"}
	Button rescaleALL,fColor=(65280,43520,0)
	Button Emptier6,pos={77,486},size={76,25},proc=Vider,title="Empty"
	Button Emptier6,help={"Removes every isotopical simulation from the agregator (top above)"}
	Button Emptier6,fColor=(65280,0,0)
	PopupMenu popup0,pos={3,512},size={149,21},bodyWidth=104,proc=PopMenuAddFromMolbags,title="Add from"
	PopupMenu popup0,mode=3,popvalue="SaveMolbag",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	PopupMenu popup1,pos={3,534},size={149,21},bodyWidth=109,proc=PopMenuSavetoMolbag,title="Save to"
	PopupMenu popup1,mode=4,popvalue="_New",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")+\"_New\""
	PopupMenu popup2,pos={2,556},size={150,21},bodyWidth=134,proc=PopMenuKillMolbag,title="Kill"
	PopupMenu popup2,mode=10,popvalue="PointsCal",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	ValDisplay valdisp0,pos={1,16},size={152,14},title="Req. resolution"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"calcreqres()"
	Button addcurrentformula,pos={72,409},size={25,12},proc=ButtonSumCurrentFormula,title="\\Z13+"
	Button addcurrentformula,help={"Changes the agregator window settings to show the entire isotopical simulation. Focus when none selected restore default settings."}
	Button addcurrentformula,fColor=(65280,21760,0)
	Button subcurrentformula,pos={72,421},size={25,12},proc=ButtonSubCurrentFormula,title="\\Z13-"
	Button subcurrentformula,help={"Changes the agregator window settings to show the entire isotopical simulation. Focus when none selected restore default settings."}
	Button subcurrentformula,fColor=(65280,21760,0)
	TitleBox title0,pos={100,404},size={73,32},variable= nommol
	TitleBox PlusOrMinus,pos={1,410},size={71,21},title="Operate to all"
	Button Delete1,pos={2,579},size={75,25},proc=ButtonProc_8,title="Erase null molecules"
	Button Delete1,fColor=(65280,0,0)
	Button Delete2,pos={2,604},size={75,25},proc=ButtonProc_9,title="Graph Current"
	Button Delete2,fColor=(65280,43520,0)
	Button Delete3,pos={78,604},size={75,25},proc=ButtonProc_10,title="Edit list"
	Button Delete3,fColor=(65280,43520,0)
	ValDisplay valdisp1,pos={2,0},size={152,14},title="Items:"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"numpnts(List_formules)"
	Button Delete4,pos={78,578},size={75,25},proc=ButtonProc_13,title="Conv. 2 ROI"
	Button Delete4,fColor=(65280,43520,0)
EndMacro

Window attributeurOLD() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(180,416,427,720)
	ListBox list0,pos={3,110},size={248,57},proc=ListPROPreturntoelaborator
	ListBox list0,help={"Display the results of decomposition by increasing delta ppm"}
	ListBox list0,fSize=10,frame=2,listWave=root:list_prop,mode= 1,selRow= 0
	ListBox list0,editStyle= 1,special= {0,15,0},widths={115,49},userColumnResize= 1
	ListBox list1,pos={2,200},size={100,106}
	ListBox list1,help={"Displays the list of patterns whose mass is used to perform linear combination to match the current m/z"}
	ListBox list1,labelBack=(65535,65535,65535),fSize=10,frame=2
	ListBox list1,listWave=root:list_molperm,mode= 1,selRow= 2,special= {0,15,0}
	SetVariable setvar0,pos={3,1},size={247,24},proc=SetVarProcInco,title="\\Z12Current m/z"
	SetVariable setvar0,help={"Displays the current mass over charge used for the calculation below. It fits to the closest datapoint automatically if Find in data is ticked."}
	SetVariable setvar0,fSize=16,format="%15.13f",fStyle=1
	SetVariable setvar0,valueBackColor=(51456,44032,58880)
	SetVariable setvar0,limits={-inf,inf,0},value= inco,live= 1
	Button button2,pos={100,26},size={80,20},proc=Buttongotonextintensity,title="Next shorter"
	Button button2,help={"Sets the current m/z to the datapoint immediatly shorter than tallest datapoint with same nominal mass as current m/z"}
	Button button2,fColor=(32768,65280,0)
	Button button3,pos={1,168},size={46,20},proc=ButtonRemoveFromPerm,title="Remove"
	Button button3,help={"Removes selected pattern from list"}
	Button button3,fColor=(65280,43520,0)
	SetVariable setvar1,pos={9,92},size={143,16},bodyWidth=60,title="Displayed results"
	SetVariable setvar1,help={"Sets the length of the proposed formulas list. Please, keep it positive."}
	SetVariable setvar1,limits={0,inf,1},value= nbprop
	ListBox list2,pos={102,200},size={150,106}
	ListBox list2,help={"Sets manual constrains to force or prevent stoechiometic abundance in proposed formulas"}
	ListBox list2,fSize=10,frame=2,listWave=root:minmax,selWave=root:driveminmax
	ListBox list2,mode= 5,special= {0,15,0}
	CheckBox check0,pos={182,170},size={63,14},title="Auto Max"
	CheckBox check0,help={"Ignores maximal manual constrains"},value= 0
	Button button4,pos={180,27},size={73,20},proc=ButtonProcheETintense,title="close and tall"
	Button button4,help={"Sets the current m/z above to the tallest datapoint that have the same nominal mass as current m/z."}
	Button button4,fColor=(32768,65280,0)
	Button button5,pos={3,50},size={80,40},proc=Buttoncalculate,title="Calculate"
	Button button5,help={"Calculates combinations of patterns (set below) that could match the current m/z (above) and sort them"}
	Button button5,fColor=(65280,43520,0)
	TitleBox title1,pos={27,186},size={35,13},title="Factors",frame=0
	TitleBox title2,pos={122,186},size={90,13},title="Manual Constraints",frame=0
	Button button0,pos={84,50},size={80,40},proc=ButtonAttribute,title="Attribute"
	Button button0,help={"Opens a panel to run automatic attributions with current settings."}
	Button button0,fColor=(65280,0,0)
	CheckBox checkFindInData,pos={4,29},size={73,14},proc=CheckFindInData_proc,title="Find in data"
	CheckBox checkFindInData,help={"Fits the current m/z to the closest datapoint after it is entered."}
	CheckBox checkFindInData,variable= findindata
	Button button6,pos={48,168},size={36,20},proc=ButtonemptyPerm,title="Empty"
	Button button6,help={"Removes selected pattern from list"},fColor=(65280,0,0)
	PopupMenu popup1,pos={84,168},size={93,21},proc=PopMenuProc_2
	PopupMenu popup1,mode=66,popvalue="ClassicalOrga",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	Button button1,pos={166,50},size={80,40},proc=ButtonProc_17,title="Graphttribution"
	Button button1,help={"Opens a panel to run automatic attributions with current settings."}
	Button button1,fColor=(52224,52224,52224)
EndMacro

Window mendeleiev() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(7,757,948,939)
	Button H,pos={1,1},size={25,25},proc=elemH,title="H",fColor=(60928,60928,60928)
	Button He,pos={807,1},size={25,25},proc=elemHe,title="He"
	Button Li,pos={1,27},size={25,25},proc=elemLi,title="Li"
	Button Be,pos={27,27},size={25,25},proc=elemBe,title="Be"
	Button B,pos={677,27},size={25,25},proc=elemB,title="B"
	Button C,pos={703,28},size={25,25},proc=elemC,title="C",fColor=(4352,4352,4352)
	Button N,pos={729,27},size={25,25},proc=elemN,title="N"
	Button N,fColor=(48896,65280,57344)
	Button O,pos={755,27},size={25,25},proc=elemO,title="O"
	Button O,fColor=(65280,48896,48896)
	Button F,pos={781,27},size={25,25},proc=elemF,title="F"
	Button Ne,pos={807,27},size={25,25},proc=elemNe,title="Ne"
	Button Na,pos={1,53},size={25,25},proc=elemNa,title="Na"
	Button Mg,pos={27,53},size={25,25},proc=elemMg,title="Mg"
	Button Al,pos={677,53},size={25,25},proc=elemAl,title="Al"
	Button Si,pos={703,53},size={25,25},proc=elemSi,title="Si"
	Button P,pos={729,53},size={25,25},proc=elemP,title="P"
	Button S,pos={755,53},size={25,25},proc=elemS,title="S"
	Button S,fColor=(65280,65280,48896)
	Button Cl,pos={781,53},size={25,25},proc=elemCl,title="Cl"
	Button Ar,pos={807,53},size={25,25},proc=elemAr,title="Ar"
	Button K,pos={1,79},size={25,25},proc=elemK,title="K"
	Button Ca,pos={27,79},size={25,25},proc=elemCa,title="Ca"
	Button Sc,pos={417,79},size={25,25},proc=elemSc,title="Sc"
	Button Ti,pos={443,79},size={25,25},proc=elemTi,title="Ti"
	Button V,pos={469,79},size={25,25},proc=elemV,title="V"
	Button Cr,pos={495,79},size={25,25},proc=elemCr,title="Cr"
	Button Mn,pos={521,79},size={25,25},proc=elemMn,title="Mn"
	Button Fe,pos={547,79},size={25,25},proc=elemFe,title="Fe"
	Button Co,pos={573,79},size={25,25},proc=elemCo,title="Co"
	Button Ni,pos={599,79},size={25,25},proc=elemNi,title="Ni"
	Button Cu,pos={625,79},size={25,25},proc=elemCu,title="Cu"
	Button Zn,pos={651,79},size={25,25},proc=elemZn,title="Zn"
	Button Ga,pos={677,79},size={25,25},proc=elemGa,title="Ga"
	Button Ge,pos={703,79},size={25,25},proc=elemGe,title="Ge"
	Button As,pos={729,79},size={25,25},proc=elemAs,title="As"
	Button Se,pos={755,79},size={25,25},proc=elemSe,title="Se"
	Button Br,pos={781,79},size={25,25},proc=elemBr,title="Br"
	Button Kr,pos={807,79},size={25,25},proc=elemKr,title="Kr"
	Button Rb,pos={1,105},size={25,25},proc=elemRb,title="Rb"
	Button Sr,pos={27,105},size={25,25},proc=elemSr,title="Sr"
	Button Y,pos={417,105},size={25,25},proc=elemY,title="Y"
	Button Zr,pos={443,105},size={25,25},proc=elemZr,title="Zr"
	Button Nb,pos={469,105},size={25,25},proc=elemNb,title="Nb"
	Button Mo,pos={495,105},size={25,25},proc=elemMo,title="Mo"
	Button Tc,pos={521,105},size={25,25},proc=elemTc,title="Tc"
	Button Ru,pos={547,105},size={25,25},proc=elemRu,title="Ru"
	Button Rh,pos={573,105},size={25,25},proc=elemRh,title="Rh"
	Button Pd,pos={599,105},size={25,25},proc=elemPd,title="Pd"
	Button Ag,pos={625,105},size={25,25},proc=elemAg,title="Ag"
	Button Cd,pos={651,105},size={25,25},proc=elemCd,title="Cd"
	Button In,pos={677,105},size={25,25},proc=elemIn,title="In"
	Button Sn,pos={703,105},size={25,25},proc=elemSn,title="Sn"
	Button Sb,pos={729,105},size={25,25},proc=elemSb,title="Sb"
	Button Te,pos={755,105},size={25,25},proc=elemTe,title="Te"
	Button I,pos={781,105},size={25,25},proc=elemI,title="I"
	Button Xe,pos={807,105},size={25,25},proc=elemXe,title="Xe"
	Button Cs,pos={1,131},size={25,25},proc=elemCs,title="Cs"
	Button Ba,pos={27,131},size={25,25},proc=elemBa,title="Ba"
	Button La,pos={53,131},size={25,25},proc=elemLa,title="La"
	Button Ce,pos={79,131},size={25,25},proc=elemCe,title="Ce"
	Button Pr,pos={105,131},size={25,25},proc=elemPr,title="Pr"
	Button Nd,pos={131,131},size={25,25},proc=elemNd,title="Nd"
	Button Pm,pos={157,131},size={25,25},proc=elemPm,title="Pm"
	Button Sm,pos={183,131},size={25,25},proc=elemSm,title="Sm"
	Button Eu,pos={209,131},size={25,25},proc=elemEu,title="Eu"
	Button Gd,pos={235,131},size={25,25},proc=elemGd,title="Gd"
	Button Tb,pos={261,131},size={25,25},proc=elemTb,title="Tb"
	Button Dy,pos={287,131},size={25,25},proc=elemDy,title="Dy"
	Button Ho,pos={313,131},size={25,25},proc=elemHo,title="Ho"
	Button Er,pos={339,131},size={25,25},proc=elemEr,title="Er"
	Button Tm,pos={365,131},size={25,25},proc=elemTm,title="Tm"
	Button Yb,pos={391,131},size={25,25},proc=elemYb,title="Yb"
	Button Lu,pos={417,131},size={25,25},proc=elemLu,title="Lu"
	Button Hf,pos={443,131},size={25,25},proc=elemHf,title="Hf"
	Button Ta,pos={469,131},size={25,25},proc=elemTa,title="Ta"
	Button W,pos={495,131},size={25,25},proc=elemW,title="W"
	Button Re,pos={521,131},size={25,25},proc=elemRe,title="Re"
	Button Os,pos={547,131},size={25,25},proc=elemOs,title="Os"
	Button Ir,pos={573,131},size={25,25},proc=elemIr,title="Ir"
	Button Pt,pos={599,131},size={25,25},proc=elemPt,title="Pt"
	Button Au,pos={625,131},size={25,25},proc=elemAu,title="Au"
	Button Hg,pos={651,131},size={25,25},proc=elemHg,title="Hg"
	Button Tl,pos={677,131},size={25,25},proc=elemTl,title="Tl"
	Button Pb,pos={703,131},size={25,25},proc=elemPb,title="Pb"
	Button Bi,pos={729,131},size={25,25},proc=elemBi,title="Bi"
	Button Po,pos={755,131},size={25,25},proc=elemPo,title="Po"
	Button At,pos={781,131},size={25,25},proc=elemAt,title="At"
	Button Rn,pos={807,131},size={25,25},proc=elemRn,title="Rn"
	Button Fr,pos={1,157},size={25,25},proc=elemFr,title="Fr"
	Button Ra,pos={27,157},size={25,25},proc=elemRa,title="Ra"
	Button Ac,pos={53,157},size={25,25},proc=elemAc,title="Ac"
	Button Th,pos={79,157},size={25,25},proc=elemTh,title="Th"
	Button Pa,pos={105,157},size={25,25},proc=elemPa,title="Pa"
	Button U,pos={131,157},size={25,25},proc=elemU,title="U"
	Button Np,pos={157,157},size={25,25},proc=elemNp,title="Np"
	Button Pu,pos={183,157},size={25,25},proc=elemPu,title="Pu"
	Button Am,pos={209,157},size={25,25},proc=elemAm,title="Am"
	Button Cm,pos={235,157},size={25,25},proc=elemCm,title="Cm"
	Button Bk,pos={261,157},size={25,25},proc=elemBk,title="Bk"
	Button Cf,pos={287,157},size={25,25},proc=elemCf,title="Cf"
	Button Es,pos={313,157},size={25,25},proc=elemEs,title="Es"
	Button Fm,pos={339,157},size={25,25},proc=elemFm,title="Fm"
	Button Md,pos={365,157},size={25,25},proc=elemMd,title="Md"
	Button No,pos={391,157},size={25,25},proc=elemNo,title="No"
	Button Lr,pos={417,157},size={25,25},proc=elemLr,title="Lr"
	Button Rf,pos={443,157},size={25,25},proc=elemRf,title="Rf"
	Button Db,pos={469,157},size={25,25},proc=elemDb,title="Db"
	Button Sg,pos={495,157},size={25,25},proc=elemSg,title="Sg"
	Button Bh,pos={521,157},size={25,25},proc=elemBh,title="Bh"
	Button Hs,pos={547,157},size={25,25},proc=elemHs,title="Hs"
	Button Mt,pos={573,157},size={25,25},proc=elemMt,title="\\M\\S14\\MN"
	Button Uun,pos={599,157},size={25,25},proc=elemUun,title="\\M\\S15\\MN"
	Button Uuu,pos={625,157},size={25,25},proc=elemUuu,title="\\M\\S1\\MH"
	Button Uub,pos={651,157},size={25,25},proc=elemUub,title="\\M\\S2\\MH"
	Button Uuq,pos={677,157},size={25,25},proc=elemUuq,title="\\M\\S12\\MC"
	Button Uuh,pos={703,157},size={25,25},proc=elemUuh,title="\\M\\S13\\MC"
	Button Zeromol,pos={832,53},size={102,50},proc=ButtonReset,title="Reset"
	Button Zeromol,fColor=(65280,43520,0)
	SetVariable setvarkrikri,pos={84,53},size={269,16},bodyWidth=45,proc=SetKriter,title="Isotopic combo probability to take into account"
	SetVariable setvarkrikri,limits={0,1,0.01},value= krikriter
	CheckBox dis_selector,pos={53,79},size={147,14},proc=Display_sel,title="Display relative to maximum"
	CheckBox dis_selector,variable= dismode
	Button adder,pos={833,1},size={102,51},proc=Ajouter,title="Add to agregator"
	Button adder,fColor=(32768,65280,0)
	SetVariable charge_setter,pos={60,105},size={91,16},bodyWidth=36,proc=RegleCharge,title="Ion charge"
	SetVariable charge_setter,limits={-5,5,1},value= charge
	Button toperm,pos={53,1},size={102,51},proc=ButtonAddtoPerm,title="Add as factor"
	Button toperm,fColor=(32768,65280,0)
	TitleBox title0,pos={313,1},size={9,9},labelBack=(65535,65535,65535)
	TitleBox title0,variable= nommol
	Button msmser,pos={833,157},size={102,25},proc=gotoMSMS,title="doMSMS"
	Button msmser,fColor=(57344,65280,48896)
	ListBox list0,pos={1,187},size={50,20}
	SetVariable setvar0,pos={159,35},size={175,16},proc=Mendev_Formula_Set,title="Set this formula:"
	SetVariable setvar0,value= _STR:"C11H22"
	SetVariable setvar1,pos={338,35},size={110,16},proc=Mendev_Formula_Add,title="or add this:"
	SetVariable setvar1,value= _STR:""
	SetVariable setvar2,pos={451,35},size={110,16},proc=Mendev_Formula_Add,title="or this:"
	SetVariable setvar2,value= _STR:""
EndMacro

Window slopecompass() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(249.75,441.5,474,613.25) rawyplode vs rawxplode
	AppendToGraph cauchyyplode vs cauchyxplode
	ModifyGraph lSize=1.2
	ModifyGraph rgb(rawyplode)=(0,0,65280)
	ModifyGraph zero(left)=1
	ModifyGraph noLabel(bottom)=2
	ModifyGraph axOffset(left)=-1.85714,axOffset(bottom)=-1.8
	ModifyGraph axThick(bottom)=0
	Label left "Slope compass"
EndMacro

Window slopehisto() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(1224.75,623.75,1401,739.25) fit_histoslo,histoslo
	ModifyGraph lSize=2
	ModifyGraph rgb(fit_histoslo)=(30464,30464,30464),rgb(histoslo)=(0,15872,65280)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=20
	ModifyGraph fStyle=1
	ModifyGraph axOffset(left)=-2
	ModifyGraph axThick=2
	Label left "Histogram"
	Label bottom "Point to point slope value"
	SetAxis bottom -0.00801870967741936,0.00802258064516129
EndMacro

Function ButtonProc_8(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			killmissmatch()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function ButtonProc_9(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			prepareCompareAttribution("Current")
			removeFakes()
			// Démarre l'interface GUI de compare Attrib
			interfaceCompareAttribution()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



////////Smooth operation

function convoluateur(ay,ax,by,bx)// a=ondelette
wave ay,ax,by,bx
variable n=numpnts(bx), m=numpnts(ax), k, j, sca = sum(ay)
ay/=sca
make/O/D/N=(n) roi0, roi1
roi1=0
	for(k=0;k<n;k+=1)
	roi0[k]=bx[k]
		for(j=0;j<m;j+=1)
			findlevel/Q bx, (bx[k]-ax[j])
			roi1[k]+=by[v_levelx]*ay[j]
		endfor
	endfor
end

function rescaler(wx,wy,step)
wave wx, wy
variable step
variable b1=wavemin(wx),b2=wavemax(wx), range=b2-b1, k, nb=range/step
make/O/D/N=(nb) roi1, roi0
for(k=0;k<nb;k+=1)
	findlevel/Q wx b1+k*step
	roi0[k]=wx[V_levelX]
	roi1[k]=wy[V_levelX]
endfor
end

function genwavelet(nb,step, type)
variable nb,step, type
nvar simu_reso,simu_shape,simu_asym, simu_mex
variable sca
nb=nb+(mod(nb,2)==0)
Make/O/D/N=(nb) waveletx, wavelety
waveletx=(x-floor(nb/2))*step
	if(type==1)
		wavelety=pearson4({1,0,simu_reso,simu_shape,simu_asym},waveletx)
	elseif(type==2)
		wavelety=sinc(waveletx/simu_mex)
	endif
sca=sum(wavelety)
wavelety/=sca
sca=wavemax(wavelety)
findlevels/Q wavelety sca/2
duplicate/O W_findlevels FWHMy, FWHMx
FWHMx=waveletx[FWHMy[x]]
FWHMy=sca/2
end

function convolutor()
convolve/A wavelety roi1
end

///end smooth operation

Function SetVarProc_4(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	wave cursint
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			cursint[0]=sva.dval
			majcursint(cursint[0])
			roifromcursint()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function SetVarProc_5(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	wave cursint
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			cursint[1]=sva.dval
			majcursint(cursint[0])
			roifromcursint()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function PopMenuProc_1(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
nvar functype, ech_wavelet, stepwidth
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			functype=pa.popNum
			killsituacontrol("Convolver","situa")
			switch(pa.popNum)
				case 1:
					SetVariable Situa1_width title="Width parameter",pos={50,78}, value=simu_reso,size={100,16},bodyWidth=50,proc=SetVarProc_genwavelet
					SetVariable Situa1_shape title="Shape parameter",pos={211,78}, value=simu_shape,size={100,16},bodyWidth=50,proc=SetVarProc_genwavelet
					SetVariable Situa1_asym title="Asym parameter",pos={356,78}, value=simu_asym,size={100,16},bodyWidth=50,proc=SetVarProc_genwavelet
					break
				case 2:
					SetVariable Situa2_width title="Width parameter",pos={50,78}, value=simu_mex,size={100,16},bodyWidth=50,proc=SetVarProc_genwavelet, limits={-inf,inf,0.0001}
					break
			endswitch
			genwavelet(ech_wavelet,stepwidth, functype)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function killsituacontrol(win,key)
string win,key
key="*"+key+"*"
string ctrlist=controlnamelist(win,";",key), lectrl
variable n=itemsinlist(ctrlist),k
for(k=0;k<n;k+=1)
	lectrl=stringfromlist(k,ctrlist)
	killcontrol/W=$win $lectrl
endfor
end

Window convolver() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1005,596,1472,993) /K=1
	GroupBox group0,pos={0,37},size={467,338},title="Wavelet"
	PopupMenu popup0,pos={11,56},size={144,21},proc=PopMenuProc_1,title="Wavelet type"
	PopupMenu popup0,mode=1,popvalue="Bell shape",value= #"\"Bell shape;Mexican hat\""
	SetVariable Perma_ech,pos={11,355},size={134,16},bodyWidth=50,proc=SetVarProc_genwavelet,title="Number of points"
	SetVariable Perma_ech,limits={0,inf,1},value= ech_wavelet
	SetVariable Perma_stepwidth,pos={160,355},size={104,16},bodyWidth=50,proc=SetVarProc_genwavelet,title="Step width"
	SetVariable Perma_stepwidth,limits={0,inf,0.001},value= stepwidth
	GroupBox datagroup,pos={0,0},size={467,37},title="Data to convolve"
	ValDisplay valdisp0,pos={11,20},size={168,14},title="Points in current data"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp0,value= #"numpnts(wave0stripped)"
	ValDisplay valdisp1,pos={182,19},size={63,14},title=", from"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"wavemin(wave0stripped)"
	ValDisplay valdisp2,pos={248,19},size={52,14},title="to"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp2,value= #"wavemax(wave0stripped)"
	ValDisplay valdisp3,pos={310,19},size={154,14},title="average stepwidth"
	ValDisplay valdisp3,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp3,value= #"(wavemax(wave0stripped)-wavemin(wave0stripped))/numpnts(wave0stripped)"
	Button button0,pos={0,376},size={467,20},proc=ButtonProc_11,title="Convolve and put into ROI"
	Button button0,fColor=(32768,65280,0)
	ValDisplay valdisp4,pos={303,220},size={154,14},title="average stepwidth"
	ValDisplay valdisp4,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp4,value= #"(wavemax(wave0stripped)-wavemin(wave0stripped))/numpnts(wave0stripped)"
	ValDisplay valdisp5,pos={272,355},size={188,14},title="points after convolution :"
	ValDisplay valdisp5,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp5,value= #"(wavemax(wave0stripped)-wavemin(wave0stripped))/stepwidth"
	SetVariable Situa1_width,pos={18,78},size={132,16},bodyWidth=50,proc=SetVarProc_genwavelet,title="Width parameter"
	SetVariable Situa1_width,value= simu_reso
	SetVariable Situa1_shape,pos={176,78},size={135,16},bodyWidth=50,proc=SetVarProc_genwavelet,title="Shape parameter"
	SetVariable Situa1_shape,value= simu_shape
	SetVariable Situa1_asym,pos={327,78},size={129,16},bodyWidth=50,proc=SetVarProc_genwavelet,title="Asym parameter"
	SetVariable Situa1_asym,value= simu_asym
	Display/W=(5,100,459,353)/HOST=#  wavelety vs waveletx
	AppendToGraph FWHMy,FWHMy vs FWHMx
	ModifyGraph mode(wavelety)=4,mode(FWHMy)=4,mode(FWHMy#1)=3
	ModifyGraph marker(wavelety)=19,marker(FWHMy)=9
	ModifyGraph rgb(FWHMy)=(0,0,65280),rgb(FWHMy#1)=(0,0,0)
	ModifyGraph mrkThick(FWHMy)=2
	ModifyGraph useMrkStrokeRGB(wavelety)=1
	ModifyGraph textMarker(FWHMy#1)={FWHMx,"default",0,90,1,0.00,50.00}
	ModifyGraph zero(bottom)=1
	ErrorBars wavelety Y,wave=(,wavelety)
	ErrorBars FWHMy Y,wave=(,FWHMy)
	RenameWindow #,G0
	SetActiveSubwindow ##
EndMacro

function initiateconvolutor()
variable/g ech_wavelet, stepwidth, functype
variable/g simu_mex
end

Function SetVarProc_genwavelet(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
nvar ech_wavelet,stepwidth, functype
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			genwavelet(ech_wavelet,stepwidth, functype)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function ButtonProc_11(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar ech_wavelet, stepwidth, functype
wave wave0stripped, wave1stripped
	switch( ba.eventCode )
		case 2: // mouse up
			rescaler(wave0stripped, wave1stripped, stepwidth)// click code here
			convolutor()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_12(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			initiateconvolutor()// click code here
			execute "convolver()"
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

////headup tree

function generatetree(nomdebag,prof,m)
variable prof,m
string nomdebag
string Smolbag="Molbag_"+nomdebag
string Schargebag="Chargebag_"+nomdebag
string elem
wave molbag=$Smolbag
wave chargebag=$Schargebag
wave molecule, mendev_masses, wave0stripped, wave1stripped,simu_mass
wave/T elements
variable n=dimsize(molbag,1),k,j, cote=(2*prof+1), t
nvar charge
Make/O/D/N=(n) lesmasses
Duplicate/O molecule transmolecule
for(k=0;k<n;k+=1)
	transmolecule=molbag[x][k]*mendev_masses[x][0]/(max(1,abs(chargebag[x])))
	lesmasses[k]=sum(transmolecule)
endfor
Make/O/D/N=(cote^n,n) cribletree
cribletree=(mod(floor(x/(cote^y)),cote)-prof)*lesmasses[y]
MatrixOP/O tree=m+sumrows(cribletree)
Duplicate/O tree, defmasstree,disttree,treeindata,treeerror,treeintens
defmasstree=defmass(tree)
disttree=tree-m
treeindata=wave0stripped[prochepic(tree,0)]
treeerror=abs(1e6*(tree-treeindata)/tree)
treeintens=wave1stripped[prochepic(tree,0)]
Make/O/D/N=(114,cote^n) Molbag_Tree=molecule[x]
for(k=0;k<n;k+=1)
	Molbag_Tree+=(mod(floor(y/(cote^k)),cote)-prof)*molbag[x][k]
endfor
Make/O/D/N=(cote^n,5) Ciblebag_Tree=0
Ciblebag_Tree[][1]=treeindata[x]
Ciblebag_Tree[][2]=treeintens[x]
Ciblebag_Tree[][3]=1e6*(tree[x]-treeindata[x]+(simu_mass[numpnts(simu_mass)-1]-m))/(tree[x]+(simu_mass[numpnts(simu_mass)-1]-m))
Make/O/D/N=(cote^n) Chargebag_Tree=charge, treeColor
Make/O/D/T/N=(cote^n) Listbag_Tree=""
for(j=0;j<cote^n;j+=1)
	for(k=0;k<114;k+=1)//verifie la presence de chaque element
			if (Molbag_Tree[k][j]!=0) // si l' element k est present
				elem=elements[k]
				Listbag_Tree[j]+=elem+"\B"+num2str(Molbag_Tree[k][j])+"\M"//genere la formule brute
			endif
		endfor
		if (charge>0)
			Listbag_Tree[j]+="\S"+num2str(abs(charge))+"+\M"
		elseif (charge<0)
			Listbag_Tree[j]+="\S"+num2str(abs(charge))+"-\M"
		endif
endfor
treecolor=(mod(floor(x/(cote^0)),cote)-prof)
Killwaves transmolecule, cribletree
end

//fitallpeaks procedure

function fittickedpeaks()
wave selpklist, peakset_current
wave/t peaklist
variable n= dimsize(selpklist,0),k
for(k=n-1;k>-1;k-=1)
	if(selpklist[k][6]==48)
		FitThisPeak(peakset_current,str2num(peaklist[k][1]),hold2cons())
		print k, peaklist[k][1]
	endif
endfor
SetTraceMixPeaks(PeakSet_Current)
end
end


function cutspec(b1,b2,guide,index,xx,yy) //coupe en incluant
wave xx,yy, guide,index
variable b1,b2
if(b1==b2)
Make/O/D/N=0 outx, outy
return -1
endif
variable bmin=min(b1,b2), bmax=max(b1,b2)
findlevel/q guide bmin
	if(numtype(V_levelx)==0)
		bmin=ceil(V_levelX)*(V_rising==1)+floor(V_levelX)*(V_rising==0)
	else
		bmin=(1-V_flag)+(V_rising==0)*numpnts(xx)
	endif
findlevel/q guide bmax
	if(numtype(V_levelx)==0)
		bmax=floor(V_levelX)*(V_rising==1)+ceil(V_levelX)*(V_rising==0)
	else
		bmax=(V_flag==1)*((V_rising==1)*numpnts(xx))
	endif
//print v_flag,v_rising,v_levelx,bmax, bmin
b1=min(bmin,bmax)
b2=max(bmin,bmax)
Make/O/D/N=(b2-b1) outx, outy
outx=xx[index[b1+x]]
outy=yy[index[b1+x]]
end

function cropbetween(b1,b2,xx,yy)//coupe un couple de waves en incluant les points entre b1 et b2 selon yy
variable b1,b2
wave xx,yy
variable bmin=min(b1,b2), bmax=max(b1,b2)
duplicate/o xx match
match=x*(yy<(bmax) && yy>(bmin))/((yy<(bmax) && yy>(bmin)))
WaveTransform zapnans match
duplicate/o match outy, outx
outy=yy[match[x]]
outx=xx[match[x]]
killwaves match
end

function calcreqres()
meltattrib("trans","")
wave trans_1, trans_0
variable n=numpnts(trans_0),k
Make/O/D/N=(n-1) resreq
resreq=abs(0.5*(trans_0[x]+trans_0[x+1])/(trans_0[x]-trans_0[x+1]))
k=wavemax(resreq)
//killwaves trans_1, trans_0
return k
end

Function ButtonProc_13(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			meltattrib("trans","")
			duplicate/O trans_1, roi1, roi2
			duplicate/O trans_0, roi0
			sort roi0, roi1, roi0
			roi2=defmass(roi0)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function tickall()
wave selpklist
selpklist[][6]=48
end

Function tickallbutton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			tickall()
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function fittickedbutton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			fittickedpeaks()
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function massloadCosimaTABS()
 	string nomfichier, fichiery, fichierx, fichierinfo
    	GetFileFolderInfo/Q
   	 NewPath/O/Q ledossier ParseFilePath(1, S_Path, ":", 1, 0)
                String list = IndexedFile(ledossier, -1, ".TAB")               
                // Sort using combined alpha and numeric sort
             //   list = SortList(list, ";", 16) 
                // Process the list
                Variable numItems = ItemsInList(list)
                Variable i
                for(i=0; i<numItems; i+=1)
                               String fileName = StringFromList(i, list)
                               Print i, fileName
                               nomfichier=fileName
   LoadWave/A=CosimaHeader/O/J/K=2/Q/D/L={0,0,93,0,0}/P=ledossier nomfichier   
    LoadWave/A=CosimaData/O/G/Q/D/L={0,0,0,1,2}/P=ledossier nomfichier   
    wave/T cosimaheader0
    wave cosimadata0, cosimadata1
    cosimaheader0=replacestring(" ",cosimaheader0,"")
    cosimaheader0=replacestring("\"",cosimaheader0,"")
    string target="",coordx="",coordy="", name="", prodid=""
    variable n=numpnts(cosimaheader0),k
    for(k=0;k<n;k+=1)
        target+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_ID",cosimaheader0[k], "=")
        coordx+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_X",cosimaheader0[k], "=")
        coordy+=stringbykey( "ROSETTA:COSIMA_SUBSTRATE_Y",cosimaheader0[k], "=")
        prodid+=stringbykey( "PRODUCT_ID",cosimaheader0[k], "=")
    endfor
    name=prodid+"_"+target+"_"+coordx+"_"+coordy
    name=name[3,22]
    string loadoption
    prompt name, "Name this Sample"
    prompt loadoption, "On-the-fly loading options", popup, "Load it as it is;Delete 0s intensities;Keep >0.5 masses"
    doprompt "Enter a good name for this sample",name, loadoption
    if(v_flag==0)
        strswitch(loadoption)
            case "Load it as it is":
               
                break
            case "Delete 0s intensities":
                duplicate/O cosimadata0 cropcosi
                cropcosi=cosimadata0/cosimadata0*x
                wavetransform zapnans cropcosi
                duplicate/O cropcosi, cropcosi0, cropcosi1
                cropcosi0=cosimadata0[cropcosi[x]]
                cropcosi1=cosimadata1[cropcosi[x]]
                duplicate/O cropcosi0 cosimadata0
                duplicate/O cropcosi1 cosimadata1
                killwaves cropcosi, cropcosi0, cropcosi1
                nomfichier+="z"
                break
            case "Crop only >1 intensities":
                break
            case "Keep >0.5 masses" :
                deletepoints 0,1, cosimadata0, cosimadata1   
                findlevel/Q cosimadata1, 0.5
                deletepoints 0,V_levelx, cosimadata0, cosimadata1
                cosimadata0=max(cosimadata0,0)
                break
        endswitch
    nomfichier=name[0,19]
    fichierinfo=nomfichier+"_info"
    fichiery=nomfichier+"_1stripped"
    fichierx=nomfichier+"_0stripped"
    Rename cosimadata1, $fichierx
    Rename cosimadata0, $fichiery
    Rename cosimaheader0, $fichierinfo
    else
    killwaves cosimaheader0, cosimadata1, cosimadata0
    endif
                endfor
                listelesdonnees()
End

function derivedpeakparam()
string lespeaksY=wavelist("Peak_*_Y",";",""), lespeaksX=wavelist("Peak_*_X",";","")
string lenomy, lenomx
variable n=itemsinlist(lespeaksY), k
make/O/D/N=(n) peakresos, peakmasses, peakFWHM
for(k=0;k<n;k+=1)
lenomy=stringfromlist(k,lespeaksY)
lenomx=stringfromlist(k,lespeaksX)
wave transy=$lenomy, transx=$lenomx
findlevels/Q transy wavemax(transy)/2
wave w_findlevels
duplicate/O w_findlevels FWHMmarkY, FWHMmarkX
FWHMmarkY=wavemax(transy)/2
FWHMmarkX=transx[w_findlevels[x]]
peakFWHM[k]=abs(FWHMmarkX[0]-FWHMmarkX[1])
peakresos[k]=mean(FWHMmarkX)/abs(FWHMmarkX[0]-FWHMmarkX[1])
peakmasses[k]=mean(FWHMmarkX)
print abs(FWHMmarkX[0]-FWHMmarkX[1]),mean(FWHMmarkX)/abs(FWHMmarkX[0]-FWHMmarkX[1])
endfor
end

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Djevahirdjian Léo
// été 2014
//
//		Transformée Hough
//			Interface de représentation de transformée de Hough
//	(assimilée transformée de Hough) pour une étude des relations
//	entre les points dans le jeu de données.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Transformée Hough
//
//		I. Parties Obligatoires a intégrer dans Attributor
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// etudeTransfoHough()
//	Initialisation de l'interface, mise en place des variables globales.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function etudeTransfoHough()
	Variable /G intervalleHoughSlopes = 0.02
	Variable /G borneHoughSlopes1 = -0.01
	Variable /G borneHoughSlopes2 = 0.01
	Variable /G intervalleHoughIntercepts = 12
	Variable /G borneHoughIntercepts1 = -6
	Variable /G borneHoughIntercepts2 = 6
	Variable /G storageHoughResolution = 600
	
	Variable /G storageBoolHoughCouleur = 0
	Variable /G storageBoolHoughReverse = 0
	Variable /G storageBoolHoughLog= 0
	
	Variable /G memoireSelectHoughTop = -1
	Variable /G memoireSelectHoughBottom = -1
	Variable /G memoireSelectHoughLeft = -1
	Variable /G memoireSelectHoughRight = -1
	
	Variable /G resoEnMasse = (intervalleHoughSlopes)/storageHoughResolution
	Variable /G posXSlope = -1
	Variable /G posYIntercept = -1
	Variable /G precslo=abs(borneHoughSlopes1-borneHoughSlopes2)/4,precinter=abs(borneHoughIntercepts1-borneHoughIntercepts2)/8
	
	Wave perm
	Wave/T list_molperm
	Duplicate /O /D perm listeMasseFacteursRosePentes
	Duplicate /O /T list_molperm listeFacteursRosePentes

	affichageTransfoHough()
	ModifyGraph/W=dmvm hideTrace(currentseg)=0
End //etudeTransfoHough()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// affichageTransfoHough()
//	Creation de l'interface d'étude de la transformée de Hough
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function affichageTransfoHough()
	
	NVar storageHoughResolution
	
	NVar storageBoolHoughCouleur
	NVar storageBoolHoughReverse
	NVar storageBoolHoughLog
	
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	NVar resoEnMasse
	NVar posXSlope
	NVar posYIntercept
	
	If(winType("etudeTransformationHough")!=0)
		KillWindow etudeTransformationHough
	endIf
	If(winType("zoomht")!=0)
		KillWindow zoomht
	endIf
	
	NewPanel/K=1 /N=etudeTransformationHough /W=(10,10,1260,720)
	Movewindow/W=etudeTransformationHough 5.25,  42.5,  939.75,  539.75
	SetWindow etudeTransformationHough,hook(imageResultatHough)=hookImageResultatHough
	// Interface réglage parametre transformée Hough (résolution"
	GroupBox conteneurParaHough win=etudeTransformationHough,title="Hough parameters", size={130,100}, pos={710,5},fstyle=1
	TitleBox titreResolutioHough win=etudeTransformationHough,title="Hough resolution", size={120,20}, pos={715,25},frame=0
	SetVariable setResolution win=etudeTransformationHough,title=" ", size={120,20}, pos={715,45},variable=storageHoughResolution,limits={50,10000,50},proc=changeHTreso
	TitleBox titreResolutioMassHough win=etudeTransformationHough,title="Mass resolution", size={120,20}, pos={715,65},frame=0
	ValDisplay affResoEnMasse win=etudeTransformationHough,title=" ",variable=resoEnMasse, size={120,20}, pos={715 ,85},bodyWidth=120
	// Interface des zooms
	GroupBox conteneurZOOM win=etudeTransformationHough,title="Tools ZOOM", size={130,100}, pos={710,115},fstyle=1
	Button boutonZOOM win=etudeTransformationHough,title="ZOOM Selection", size={120,25}, pos={715,135}, fColor=(20000,45000,3000),proc=zoomSelection
	Button boutonDEZOOM win=etudeTransformationHough,title="deZOOM Selection", size={120,25}, pos={715,160}, fColor=(65000,45000,3000),proc=dezoomSelection
	Button boutonREINIT win=etudeTransformationHough,title="ZOOM ReInitialisation", size={120,25}, pos={715,185}, fColor=(40000,20000,30000),proc=zoomReinit
	// Interface des visualisations, options sur la tansformée de Hough et apel d'autres représentations (rose de pentes)
	GroupBox conteneurVisu win=etudeTransformationHough,title="Vizualisation", size={130,140}, pos={710,225},fstyle=1
	CheckBox checkCouleur win=etudeTransformationHough,title="Color",variable=storageBoolHoughCouleur, size={100,15}, pos={715,247},proc=changeStyleGraph
	CheckBox checkReverse win=etudeTransformationHough,title="Inverted colors",variable=storageBoolHoughReverse, size={100,15}, pos={715,265},proc=changeStyleGraph
	CheckBox checkLog win=etudeTransformationHough,title="Log Scale",variable=storageBoolHoughLog, size={100,15}, pos={715,283},proc=changeStyleGraph
	Button boutonRosePentes win=etudeTransformationHough,title="Rose des Pentes", size={120,30}, pos={715,300}, fColor=(20000,45000,23000),proc=tracerRosedesPentes
	// Choix des facteurs apparaissant sur la rose des pentes ; choix d'un molbag
	PopupMenu choixMolBag win=etudeTransformationHough,title=" ", value=replacestring("Molbag_",wavelist("Molbag_*",";",""),""),pos={715,335},size={120,25},bodyWidth=120,proc=chargerMolbagFacteurAff
	// Interface affichage de la position du curseur dans la transformée de Hough
	GroupBox conteneurPos win=etudeTransformationHough,title="Position curseur", size={400,40}, pos={845,5},fstyle=1
	ValDisplay affPosYIntercept win=etudeTransformationHough,title="Intercept ",variable=posYIntercept, size={180,20}, pos={850,25},bodyWidth=90,fstyle=1
	ValDisplay affPosXSlope win=etudeTransformationHough,title="Slope ",variable=posXSlope, size={180,20}, pos={1040,25},bodyWidth=90,fstyle=1
	
	Make /O /N=(storageHoughResolution,storageHoughResolution) resultatHough2=0
	Make/O/D/N=1 intermarks, slomarks
	
	// Affichage de l'image de transformée de Hough contenue dans la wave resultatHough2
	Display/HOST=etudeTransformationHough /N=imageResultatHough /W=(5,5,705,705)
	AppendImage /W=etudeTransformationHough#imageResultatHough resultatHough2
	AppendToGraph /W=etudeTransformationHough#imageResultatHough plotinterdmvm vs plotslodmvm
	AppendToGraph /W=etudeTransformationHough#imageResultatHough intermarks vs slomarks
	ModifyGraph lstyle(intermarks)=1
	ModifyGraph/W=etudeTransformationHough#imageResultatHough mode(plotinterdmvm)=3,marker(plotinterdmvm)=55
	ModifyGraph/W=etudeTransformationHough#imageResultatHough rgb(plotinterdmvm)=(0,52224,52224),useMrkStrokeRGB(plotinterdmvm)=1
	
	// Génération des Labels des axes adapté à l'affichage pour les différents graphs
	Make /O /D /N=7 conteneurTicksIntercepts
	Make /O /T /N=7 conteneurTicksLabelsIntercepts
	Make /O /D /N=7 conteneurTicksSlopes
	Make /O /T /N=7 conteneurTicksLabelsSlopes
	Make /O /D /N=7 conteneurTicksInterceptsHisto
	Make /O /T /N=7 conteneurTicksLabelsInterHisto
	Make /O /D /N=7 conteneurTicksSlopesHisto
	Make /O /T /N=7 conteneurTicksLabelsSlopesHisto
	Wave conteneurTicksIntercepts
	Wave /T conteneurTicksLabelsIntercepts
	Wave conteneurTicksSlopes
	Wave /T conteneurTicksLabelsSlopes
	Wave conteneurTicksInterceptsHisto
	Wave /T conteneurTicksLabelsInterHisto
	Wave conteneurTicksSlopesHisto
	Wave /T conteneurTicksLabelsSlopesHisto
	Variable j = 0
	For(j=0 ; j<7 ; j+=1)
		conteneurTicksIntercepts[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsIntercepts[j] = num2str(borneHoughIntercepts1 + conteneurTicksIntercepts[j] * (borneHoughIntercepts2-borneHoughIntercepts1)/storageHoughResolution)[0,7]
		
		conteneurTicksSlopes[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsSlopes[j] = num2str(borneHoughSlopes1 + conteneurTicksSlopes[j] * (borneHoughSlopes2-borneHoughSlopes1)/storageHoughResolution)[0,7]
		
		conteneurTicksInterceptsHisto[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsInterHisto[j] = num2str(borneHoughIntercepts1 + conteneurTicksInterceptsHisto[j] * (borneHoughIntercepts2-borneHoughIntercepts1)/storageHoughResolution)[0,7]
		
		conteneurTicksSlopesHisto[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsSlopesHisto[j] = num2str(borneHoughSlopes1 + conteneurTicksSlopesHisto[j] * (borneHoughSlopes2-borneHoughSlopes1)/storageHoughResolution)[0,7]
	endFor
	
	ModifyGraph margin(left)=80,margin(bottom)=80,margin(top)=20,margin(right)=20, userticks(left)={conteneurTicksIntercepts,conteneurTicksLabelsIntercepts},fsize=10, userticks(bottom)={conteneurTicksSlopes,conteneurTicksLabelsSlopes}
	Label left "\Z14Intercepts (SI)"
	Label bottom "\Z14Slopes (SI)"
	
	// Affichage de la croix curseur vert en supperposition sur l'image
	Make /O /N=(storageHoughResolution) histoHoughIntercepts = 0
	Make /O /N=(storageHoughResolution) histoHoughSlopes = 0
	Make /O /N=2 histoHoughAffichageVert1 = storageHoughResolution
	Make /O /N=2 histoHoughAffichageVert2 = {0,storageHoughResolution}
	Make /O /N=2 histoHoughAffichageHoriz1 = storageHoughResolution
	Make /O /N=2 histoHoughAffichageHoriz2 = {0,storageHoughResolution}
	AppendToGraph /W=etudeTransformationHough#imageResultatHough /C=(15000,50000,15000) histoHoughAffichageVert2 vs histoHoughAffichageVert1
	AppendToGraph /W=etudeTransformationHough#imageResultatHough /C=(15000,50000,15000) histoHoughAffichageHoriz1 vs histoHoughAffichageHoriz2
	
	// Affichage des histos
	Display/HOST=etudeTransformationHough /N=histoInterceps /W=(845,49,1245,377) histoHoughIntercepts
	ModifyGraph /W=etudeTransformationHough#histoInterceps mode[0]=5, userticks(bottom)={conteneurTicksInterceptsHisto,conteneurTicksLabelsInterHisto},margin(left)=30,notation(left)=1, btLen(left)=2,margin(right)=5
	TextBox /A=MT /W=etudeTransformationHough#histoInterceps /F=0 /B=1 "Histo Intercepts (vertical)"
	Display/HOST=etudeTransformationHough /N=histoSlopes /W=(845,377,1245,705) histoHoughSlopes
	ModifyGraph /W=etudeTransformationHough#histoSlopes mode[0]=5, userticks(bottom)={conteneurTicksSlopesHisto,conteneurTicksLabelsSlopesHisto},margin(left)=30,notation(left)=1, btLen(left)=2,margin(right)=5
	TextBox /A=MT /W=etudeTransformationHough#histoSlopes /F=0 /B=1 "Histo Slopes (horizontal)"
	
	Display/N=zoomht /K=1/W=(690.75,335.75,1152,537.5)/T roiHTinter vs roiHTslo
	AppendImage/T resultatHough2
	ModifyImage resultatHough2 ctab= {*,*,Grays,0}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph mirror=2
	ModifyGraph nticks=12
	ModifyGraph minor=1
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2

	// Calcul transfo de Hough
	miseAjourTransfo()
	
	// Enregistre la première transformée de Hough pour génération de la rose des pentes complète
	Duplicate /O resultatHough2 resultatHough2Initial
End // affichageTransfoHough()

////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Fonction de gestion des événements de la fenêtre de transfo de
//	Hough, est apellée à chaque fosi qu'il se passe quelquechose dans
//	la fen^tre.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function hookImageResultatHough(s)
	STRUCT WMWinHookStruct &s
	
	NVar intervalleHoughSlopes
	NVar intervalleHoughIntercepts
	
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	NVar memoireSelectHoughTop
	NVar memoireSelectHoughBottom
	NVar memoireSelectHoughLeft
	NVar memoireSelectHoughRight
	
	NVar storageHoughResolution
	NVar posXSlope
	NVar posYIntercept
	nvar resoenmasse,inco,precinter,precslo
	
	Wave currentseg

	Variable hookResult = 0
	
	// On gère la fermeture de la fenetre
	If(s.eventCode == 2)
		// Fermeture fenetre
		If(Wintype("RoseDesPentes")!=0)
			KillWindow RoseDesPentes
		endIf
		If(winType("zoomht")!=0)
		KillWindow zoomht
		endIf
		//ModifyGraph/W=dmvm hideTrace(currentseg)=1
		nettoyageHough()
		return 0
	endIf
	
	// Verif si la souris est au bon endroit
	// La suite de la fonction s'executera si et seulement si la souris se trouve dans la zone de l'image 600px par 600px
	// En cas de redimensionnement de la fenêtre le fonctionnement ne peut être garanti car je n'ai pas eu le temps de coder
	// la taille de l'image de façon relative à la taille de la fenêtre (avec un hook de redimensionnement comme pour l'interface
	// Graphs de compare attrib.
	If(abs(Axisvalfrompixel("etudeTransformationHough#imageResultatHough","bottom",s.mouseLoc.h)-storageHoughResolution/2)>storageHoughResolution/2 || abs(Axisvalfrompixel("etudeTransformationHough#imageResultatHough","left",s.mouseLoc.v)-storageHoughResolution/2)>storageHoughResolution/2)
		// On est en dehors de la zone on ne fait rien.
		return 0
	endIf
	GetMarquee /W=etudeTransformationHough#imageResultatHough left,bottom
	switch(s.eventCode)
		case 3:
			//Click !
			// On récupère l'éventuelle sélection sans la tuer car on a besoin pour la suite
			If(V_Flag==0)
				print s.eventmod,hookresult
				// pas de selection
				If(s.eventMod==16 | s.eventMod==17)// c'est un click droit
					PopupContextualMenu "ZOOM;ReinitZOOM;ROI from red parallelogram;Attribute the major segment from red parallelogram"
					Switch(V_flag)
						case -1:
							// pas de selection dans menu
						Break
						case 1:
							// zoom dans selection
							zoomHTinmac()
						break
						case 2:
							// reinit zoom
							zoomReinit("")
						break
						case 3:
							// met les points qui satisfont les paramètres de slection dans la roi
							traceroiXpress(posYintercept,posXslope,precinter,precslo,inco)
						break
						case 4:
							// met les points qui satisfont les paramètres de slection dans la roi et attribue le segment le plus frequent
							gethoughmatch()
							findmostzero()
							analysemult(1)
						break
					endSwitch
					hookResult = 1
				elseif(s.eventmod==3)
					zoomHTinmac()
				endIf
			Else			
				If(s.eventMod==16  | s.eventMod==17)
					// c'est un click droit
					hookResult = 1
					PopupContextualMenu "ZOOM;deZOOM;ReinitZOOM;ROI from red parallelogram;Attribute the major segment from red parallelogram"
					Switch(V_flag)
						case -1:
							// pas de selection dans menu
						Break
						case 1:
							// zoom dans selection
							zoomSelection("")
						break
						case 2:
							// dezoom dans selection
							dezoomSelection("")
						break
						case 3:
							// reinit zoom
							zoomReinit("")
						break
						case 4:
							// met les points qui satisfont les paramètres de slection dans la roi
							gethoughmatch()
						break
						case 5:
							// met les points qui satisfont les paramètres de slection dans la roi et attribue le segment le plus frequent
							gethoughmatch()
							findmostzero()
							analysemult(1)
						break
					endSwitch
				endIf
				memoireSelectHoughTop = -1
				memoireSelectHoughBottom = -1
				memoireSelectHoughLeft = -1
				memoireSelectHoughRight = -1
			endIf
			
			break
		case 4:
			//Souris bouge -> il faut changer les histogrammes
			// On ne laisse pas Igor traiter le hook, on pourrait cela ne poserait normalement pas de problème
			//depuis igor 8 on laisse igor traiter le hook
			//hookResult =1
			if(v_flag)
				posYintercept=pix2inter((v_bottom+v_top)/2)
				posXslope=pix2slo((v_right+v_left)/2)
				precslo=(pix2slo(v_right)-pix2slo(v_left))/2
			else
			Variable indexSlopesSurvol = Axisvalfrompixel("etudeTransformationHough#imageResultatHough","bottom",s.mouseLoc.h)
			Variable indexInterceptsSurvol = Axisvalfrompixel("etudeTransformationHough#imageResultatHough","left",s.mouseLoc.v)
			Wave SlopeCursorY, SlopeCursorX, cauchyxplode,cauchyyplode		
			posXSlope = borneHoughSlopes1 + indexSlopesSurvol * (intervalleHoughSlopes/storageHoughResolution)
			posYIntercept = borneHoughIntercepts1 + indexInterceptsSurvol * (intervalleHoughIntercepts/storageHoughResolution)
			endif
			miseAJourHistoRose(indexSlopesSurvol, indexInterceptsSurvol)
			getaxis/W=dmvm/Q bottom
			currentseg[][1]={v_min, v_max}
			currentseg[][0]={posXSlope*v_min+posYIntercept,posXSlope*v_max+posYIntercept}
				If(Wintype("RoseDesPentes")!=0)
					// Si rose des pentes on la met à jour
					Wave SlopeCursorY, SlopeCursorX, rosepentesX,rosepentesY
					FindLevel /P /Q cauchyx, posxslope
					SlopeCursorX[1]=rosepentesX[V_LevelX]
					SlopeCursorX[0]=wavemin(wave0stripped)
					SlopeCursorY[1]=rosepentesY[V_LevelX]
			endif
			Break
		case 22 :
			// la molette est utilisée
			// =1 pour eviter toutes autres actions d'Igor
			hookResult =1
			//print s.eventmod
			If(s.wheelDy > 0)
				if(s.eventmod==8)
					precslo*=1.1
				elseif(s.eventmod==0)
					precinter*=1.1
				endif
			Elseif(s.wheelDy < 0)
				if(s.eventmod==8)
					precslo*=0.9
				elseif(s.eventmod==0)
					precinter*=0.9
				endif
			endIf
			hookResult =1
			break
	endswitch
	CalcHTinterGaps(posYintercept,posXslope)

	return hookResult		// 0 if nothing done, else 1
End

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Fonctions de Mise à jour des différentes fenêtres
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// miseAjourTransfo()
//	Recalcule la transformée de Hough a partir des parametres d'intervalles
//	des variables globales.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function miseAjourTransfo()
	NVar intervalleHoughSlopes
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	
	NVar intervalleHoughIntercepts
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	NVar storageHoughResolution
	Nvar precinter,precslo
	
	Wave intercepts = intercepts
	Wave slopes = slopes
	
	precslo=abs(borneHoughSlopes1-borneHoughSlopes2)/4
	precinter=abs(borneHoughIntercepts1-borneHoughIntercepts2)/8
	
	Duplicate /O intercepts intercepts3
	Duplicate /O slopes slopes3
	
	intercepts3 = round(((intercepts3 - borneHoughIntercepts1)*storageHoughResolution)/(intervalleHoughIntercepts))
	slopes3 = round(((slopes3 - borneHoughSlopes1)*storageHoughResolution)/(intervalleHoughSlopes))
	
	Make /O /N=(storageHoughResolution,storageHoughResolution) resultatHough2=0
	Wave resultatHough2
	
	Variable i = 0
	Variable j = 0
	Variable imax = dimSize(slopes3,0)
	Variable jmax = dimSize(slopes3,1)
	For(i=0 ; i< imax ; i+=1)
		For(j=0 ; j< imax ; j+=1)
			Multithread resultatHough2[slopes3[i][j]][intercepts3[i][j]] += 1
		endFor
	endFor
	
	Multithread resultatHough2[0][] = 0
	Multithread resultatHough2[][0] = 0
	Multithread resultatHough2[imax-1][] = 0
	Multithread resultatHough2[][jmax-1] = 0
	
	Wave conteneurTicksIntercepts
	Wave /T conteneurTicksLabelsIntercepts
	Wave conteneurTicksSlopes
	Wave /T conteneurTicksLabelsSlopes
	For(j=0 ; j<7 ; j+=1)
		conteneurTicksIntercepts[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsIntercepts[j] = num2str(borneHoughIntercepts1 + conteneurTicksIntercepts[j] * (borneHoughIntercepts2-borneHoughIntercepts1)/storageHoughResolution)[0,7]
		
		conteneurTicksSlopes[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsSlopes[j] = num2str(borneHoughSlopes1 + conteneurTicksSlopes[j] * (borneHoughSlopes2-borneHoughSlopes1)/storageHoughResolution)[0,7]
	endFor

	miseajourHTmarks()
	
	NVar resoEnMasse
	resoEnMasse = (intervalleHoughSlopes)/storageHoughResolution // réso en masse = réso en slope (démontré par FROD)
End //miseAjourTransfo()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// miseAJourHistoRose(indexSlopes, indexIntercepts)
//	Fonction mettant à jour les histogrammes ainsi que la rose des
//	pentes, si elle existe. Se base sur la position de la souris mais aussi
//	sur la selection en cours pour proposer des histogrammes moyennés.
//	Si la rose des pentes existe on met à jour l'affichage de la selection
//	en cours, tracé vert sur le graph.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function miseAJourHistoRose(indexSlopes, indexIntercepts)
	Variable indexSlopes
	Variable indexIntercepts
	
	NVar memoireSelectHoughTop
	NVar memoireSelectHoughBottom
	NVar memoireSelectHoughLeft
	NVar memoireSelectHoughRight
	
	// recup Selection
	GetMarquee /W=etudeTransformationHough#imageResultatHough left,bottom
	
	Wave conteneurTicksInterceptsHisto
	Wave /T conteneurTicksLabelsInterHisto
	Wave conteneurTicksSlopesHisto
	Wave /T conteneurTicksLabelsSlopesHisto
			
	NVar intervalleHoughSlopes
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	
	NVar intervalleHoughIntercepts
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	NVar storageHoughResolution

	Wave resultatHough2
	Wave histoHoughIntercepts
	Wave histoHoughSlopes
	
	// On déplace le viseur à la nouvelle position sur l'image de la transformée de Hough
	Wave histoHoughAffichageVert1
	Wave histoHoughAffichageVert2
	Wave histoHoughAffichageHoriz1
	Wave histoHoughAffichageHoriz2
	histoHoughAffichageVert1[0] = indexSlopes
	histoHoughAffichageVert1[1] = indexSlopes
	histoHoughAffichageVert2[1] = storageHoughResolution
	histoHoughAffichageHoriz1[0] = indexIntercepts
	histoHoughAffichageHoriz1[1] = indexIntercepts
	histoHoughAffichageHoriz2[1] = storageHoughResolution
	
	Variable j
	
	If(V_Flag==0)
		// Pas de selection en cours
		// ON présente un histogramme de position (on fait des tranches correspondant au viseur vert sur la transformée de hough

		// Récup' des histogrammes
		Redimension /N=(storageHoughResolution) histoHoughIntercepts, histoHoughSlopes
		//For(j=0 ; j<storageHoughResolution ; j+=1)
			Multithread histoHoughIntercepts=resultatHough2[indexSlopes][x]
			Multithread  histoHoughSlopes=resultatHough2[x][indexIntercepts]
		//endFor
		
		//Fabrique les labels de axes pour correspondre à l'actuel zoom
		For(j=0 ; j<7 ; j+=1)
		conteneurTicksInterceptsHisto[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsInterHisto[j] = num2str(borneHoughIntercepts1 + conteneurTicksInterceptsHisto[j] * (intervalleHoughIntercepts)/storageHoughResolution)[0,7]
		
		conteneurTicksSlopesHisto[j] = 0 + j * (storageHoughResolution/6)
		conteneurTicksLabelsSlopesHisto[j] = num2str(borneHoughSlopes1 + conteneurTicksSlopesHisto[j] * (intervalleHoughSlopes)/storageHoughResolution)[0,7]
		endFor
	Else
		// Selection en cours, on vérifie sir la selection a changée
		If(boolSelectionChanged()==1)
			// Selection différente de celle en mémoire -> il faut afficher histogrammes de selection		
			// Selection en cours on fait l'histogramme dans la selection
			Duplicate /O /R=[V_left,V_right][V_bottom,V_top] resultatHough2, resultatHough2Selection
			Redimension /N=(DimSize(resultatHough2Selection,1)) histoHoughIntercepts
			Redimension /N=(DimSize(resultatHough2Selection,0)) histoHoughSlopes
			Make /O /N=(DimSize(resultatHough2Selection,1)) memoireTestHough
		
			// On écrase la selection sur les axes pour générer nos fonction d'histo
			Matrixop/O memoireTestHough=sumCols(resultatHough2Selection)
			Matrixop/O histoHoughSlopes=sumRows(resultatHough2Selection)
			Variable i
			For(i=0 ; i<numpnts(memoireTestHough) ; i+=1)
				histoHoughIntercepts[i] = memoireTestHough[0][i]
			endFor
			
			// On génère les labels d'axes des histos correspondant à notre selection
			j = 0
			For(j=0 ; j<7 ; j+=1)
				conteneurTicksInterceptsHisto[j] = 0 + j * (abs(V_top-V_bottom)/6)
				conteneurTicksLabelsInterHisto[j] = num2str(borneHoughIntercepts1+( (V_bottom + conteneurTicksInterceptsHisto[j]) * intervalleHoughIntercepts/storageHoughResolution))[0,7]
		
				conteneurTicksSlopesHisto[j] = 0 + j * (abs(V_right-V_left)/6)
				conteneurTicksLabelsSlopesHisto[j] = num2str(borneHoughSlopes1 + ((V_left + conteneurTicksSlopesHisto[j]) * intervalleHoughSlopes/storageHoughResolution))[0,7]
			endFor
			
			If(Wintype("RoseDesPentes")!=0)
				// Si rose des pentes on la met à jour pour prendre en compte la nouvelle selection
				miseAJourRosePentes()
			endIf
			
			Killwaves /Z resultatHough2Selection,memoireTestHough
		endIf
	endIf
End // miseAJourHistoRose(indexSlopes, indexIntercepts)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// miseAJourRosePentes()
//	Fonction mettant à jour le tracé de la rose des pentes.
//	Met à jour nottament la courbe de selection verte et les facteurs
//	de referecne, calculé à partir du molbag sélectionné.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function miseAJourRosePentes()
	NVar storageHoughResolution
	
	Wave resultatHough2Initial
	Wave histoHoughSlopes
	Wave histoHoughIntercepts
	Wave wave0stripped
	
	// Fabrique histo de l'etat initial
	Make /O /N=(DimSize(resultatHough2Initial,0)) histoHoughSlopeRoseGen
	Matrixop/O histoHoughSlopeRoseGen=sumRows(resultatHough2Initial)
	
	
	Make /O /N=(numpnts(histoHoughSlopeRoseGen)) rosePentesX,xplode
	Make /O /N=(numpnts(histoHoughSlopeRoseGen)) rosePentesY,yplode
	
	Duplicate /O histoHoughSlopeRoseGen histoHoughRoseGenSlope, histoHoughRoseGenRayon
	histoHoughRoseGenSlope = -0.01 + (x*(0.02/storageHoughResolution))
	histoHoughRoseGenRayon = sqrt(histoHoughRoseGenSlope*histoHoughRoseGenSlope + histoHoughSlopeRoseGen*histoHoughSlopeRoseGen)
	
	rosePentesX = histoHoughRoseGenRayon * cos(atan(histoHoughRoseGenSlope))
	rosePentesY = histoHoughRoseGenRayon * sin(atan(histoHoughRoseGenSlope))
	
	// Fabrique Chauchy
	CurveFit/Q /W=2 lor, histoHoughSlopeRoseGen/D
	Wave W_coef
	Make/O/N=(numpnts(histoHoughSlopeRoseGen)) cauchyX, cauchyY, cauchyxplode, cauchyyplode
	cauchyx  =histoHoughRoseGenSlope
	Variable k
	For(k=0 ; k<numpnts(cauchyy) ; k+=1)
		cauchyy[k] = W_coef[0] + W_coef[1] / ( (k-W_coef[2])^2 + W_coef[3])
	endFor
	
	cauchyxplode=cauchyy*cos(atan(cauchyx))
	cauchyyplode=cauchyy*sin(atan(cauchyx))
	
	xplode=max(0,(histoHoughRoseGenRayon-cauchyy))*cos(atan(histoHoughRoseGenSlope))
	yplode=max(0,(histoHoughRoseGenRayon-cauchyy))*sin(atan(histoHoughRoseGenSlope))
	Variable maxx=max(wavemax(xplode), abs(wavemin(xplode)))
	Variable maxy=max(wavemax(yplode), abs(wavemin(yplode)))
	Variable contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
	xplode*=contrainte
	yplode*=contrainte
	cauchyxplode*=contrainte
	cauchyyplode*=contrainte
	rosePentesX*=contrainte
	rosePentesY*=contrainte
	Wave wave0stripped
	//xplode+=wavemin(wave0stripped)
	//cauchyxplode+=wavemin(wave0stripped)
	//rosePentesX+=wavemin(wave0stripped)
	
	GetMarquee /W=etudeTransformationHough#imageResultatHough left,bottom
	//Ajoutte rose pente selection
	If(V_Flag!=0)
		NVar borneHoughSlopes1
		NVar intervalleHoughSlopes
		NVar storageHoughResolution
	
		Make /O /N=(numpnts(histoHoughSlopes)) rosePentesXSelect
		Make /O /N=(numpnts(histoHoughSlopes)) rosePentesYSelect

		Duplicate /O histoHoughSlopes histoHoughRoseSlopeSelect histoHoughRoseGenRSelect
		Variable i
		For(i=0 ; i<numpnts(histoHoughRoseSlopeSelect) ; i +=1)
			histoHoughRoseSlopeSelect[i] = borneHoughSlopes1 + ((V_left + i) * intervalleHoughSlopes/storageHoughResolution)
			histoHoughRoseGenRSelect[i] = sqrt(histoHoughRoseSlopeSelect[i]*histoHoughRoseSlopeSelect[i] + histoHoughSlopes[i]*histoHoughSlopes[i])
		endFor
	
		rosePentesXSelect = histoHoughRoseGenRSelect * cos(atan(histoHoughRoseSlopeSelect))
		rosePentesYSelect = histoHoughRoseGenRSelect * sin(atan(histoHoughRoseSlopeSelect))
		
		rosePentesXSelect*=contrainte
		rosePentesYSelect*=contrainte
		//rosePentesXSelect+=wavemin(wave0stripped)
	endIf
	
	//Affichage des facteurs d'attribution
	Wave listeMasseFacteursRosePentes
	Wave/T listeFacteursRosePentes
	// NETTOYAGE des anciens facteurs
	Variable indexCompteur = 0
	
	Variable boolTestPresenceFonction = 1
	Do
		String nomWaveASupprimer1 =  "roseDesPentesFacteur" + num2str(indexCompteur) + "X"
		String nomWaveASupprimer2 =  "roseDesPentesFacteur" + num2str(indexCompteur) + "Y"
		boolTestPresenceFonction = waveexists($nomWaveASupprimer2)
		
		If(boolTestPresenceFonction !=0)		
			RemoveFromGraph /W=RoseDesPentes $nomWaveASupprimer2
			Tag /C /N=$nomWaveASupprimer2 /W=RoseDesPentes /K
			DoUpdate
			KillWaves /Z $nomWaveASupprimer1 $nomWaveASupprimer2
			
			indexCompteur+=1
		endIf
	While(boolTestPresenceFonction != 0)
	// Nouveaux facteurs
	Variable l = 0
	For(l=0 ; l<numpnts(listeMasseFacteursRosePentes) ; l+=1)
		String nomWaveX = "roseDesPentesFacteur" + num2str(l) + "X"
		String nomWaveY = "roseDesPentesFacteur" + num2str(l) + "Y"
		Make /O /N=2 $nomWaveX
		Make /O /N=2 $nomWaveY
		Wave roseEnConstructionX  = $nomWaveX
		Wave roseEnConstructionY  = $nomWaveY
		Variable penteEnconstruction = (listeMasseFacteursRosePentes[l] - round(listeMasseFacteursRosePentes[l])) / listeMasseFacteursRosePentes[l]
		FindLevel /P /Q cauchyx, penteEnconstruction
		roseEnConstructionX[0]=0//+wavemin(wave0stripped)
		roseEnConstructionX[1]=rosepentesX[V_LevelX]//-wavemin(wave0stripped)
		roseEnConstructionY[0]=0
		roseEnConstructionY[1]=rosepentesY[V_LevelX]
		AppendToGraph /W=RoseDesPentes /C=(0,0,0) roseEnConstructionY vs roseEnConstructionX
		
		Tag /C /N=$nomWaveY /W=RoseDesPentes $nomWaveY,cauchyxplode[V_LevelX] ,listeFacteursRosePentes[l]
	endFor
End // miseAJourRosePentes()

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Fonctions Zooms
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// zoomSelection(nomBouton)
//	Zoom dans la zone de selection -> la zone de selection devient la
//	zone de la transformée recalculée.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function zoomSelection(nomBouton) : ButtonControl
	String nomBouton
	
	NVar intervalleHoughSlopes
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar intervalleHoughIntercepts
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	NVar storageHoughResolution
	
	// Recup zone selection
	GetMarquee /W=etudeTransformationHough#imageResultatHough /K left,bottom
	
	If(V_Flag==0)
		// Pas de selection, Ozef
	Else
		//calcul nouvelles bornes
		borneHoughSlopes2 = borneHoughSlopes1 + V_right/storageHoughResolution * intervalleHoughSlopes
		borneHoughSlopes1 = borneHoughSlopes1 + V_left/storageHoughResolution * intervalleHoughSlopes
		intervalleHoughSlopes = borneHoughSlopes2-borneHoughSlopes1
		
		borneHoughIntercepts2 = borneHoughIntercepts1 + V_top/storageHoughResolution * intervalleHoughIntercepts
		borneHoughIntercepts1 = borneHoughIntercepts1 + V_bottom/storageHoughResolution * intervalleHoughIntercepts
		intervalleHoughIntercepts = borneHoughIntercepts2-borneHoughIntercepts1
		
		// recalcul transfo de Hough
		miseAjourTransfo()
	endIf
End // zoomSelection(nomBouton)

Function zoomHTinMac()
	
	NVar intervalleHoughSlopes
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar intervalleHoughIntercepts
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	NVar storageHoughResolution
	
	wave roiHTinter, roiHTslo
		borneHoughSlopes2 = pix2slo(wavemax(roiHTslo))
		borneHoughSlopes1 = pix2slo(wavemin(roiHTslo))
		intervalleHoughSlopes = borneHoughSlopes2-borneHoughSlopes1
		
		borneHoughIntercepts2 = pix2inter(wavemax(roiHTinter))
		borneHoughIntercepts1 = pix2inter(wavemin(roiHTinter))
		intervalleHoughIntercepts = borneHoughIntercepts2-borneHoughIntercepts1
		
		// recalcul transfo de Hough
		miseAjourTransfo()
End

////////////////////////////////////////////////////////////////////////////////////////////////////////
// dezoomSelection(nomBouton)
//	Dezoom en se basant sur la zone sélectionée. Le centre devient le
//	centre de la zone selectionée et la largueur de la zone selectionnee
//	pilotte directement la taille de la transformée de Hough recalculée.
//	Plus la zone selectionnée est grande plsu le dezoom est important
//	(zone recalculée plus grande).
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function dezoomSelection(nomBouton) : ButtonControl
	String nomBouton
	
	NVar intervalleHoughSlopes
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar intervalleHoughIntercepts
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	NVar storageHoughResolution
	
	// Recup de la selection
	GetMarquee /W=etudeTransformationHough#imageResultatHough /K left,bottom
	
	If(V_Flag==0)
		// Il n'y a pas de selection
	Else
		// Calcul des nouvelles bornes
		borneHoughSlopes2 = borneHoughSlopes2 + ((abs(V_right-V_left))/storageHoughResolution * intervalleHoughSlopes)
		borneHoughSlopes1 = borneHoughSlopes1 - ((abs(V_right-V_left))/storageHoughResolution * intervalleHoughSlopes)
		intervalleHoughSlopes = borneHoughSlopes2-borneHoughSlopes1
		
		borneHoughIntercepts2 = borneHoughIntercepts2 + ((abs(V_top-V_bottom))/storageHoughResolution * intervalleHoughIntercepts)
		borneHoughIntercepts1 = borneHoughIntercepts1 - ((abs(V_top-V_bottom))/storageHoughResolution * intervalleHoughIntercepts)
		intervalleHoughIntercepts = borneHoughIntercepts2-borneHoughIntercepts1
		
		// Recalcul de la transfo de Hough
		miseAjourTransfo()
	endIf
End // dezoomSelection(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// zoomReinit(nomBouton)
//	Reinitialise le zoom à la valeur initiale en réinitialisant les bornes
//	de calcul de la transformée.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function zoomReinit(nomBouton) : ButtonControl
	String nomBouton
	
	Variable /G intervalleHoughSlopes = 0.02
	Variable /G borneHoughSlopes1 = -0.01
	Variable /G borneHoughSlopes2 = 0.01
	Variable /G intervalleHoughIntercepts = 12
	Variable /G borneHoughIntercepts1 = -6
	Variable /G borneHoughIntercepts2 = 6
	
	// Recup de la selection pour l'effacer
	GetMarquee /W=etudeTransformationHough#imageResultatHough /K left,bottom
	
	miseAjourTransfo()
	miseAJourHistoRose(borneHoughSlopes1,borneHoughIntercepts1)
End

////////////////////////////////////////////////////////////////////////////////////////////////////////
// zoomPar2(posX,posY)
//	Recalcule la transformée de Hough pour avoir un zoom par deux,
//	zone centrée sur posX posY deux fois moins large
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function zoomPar2(posX,posY)
	Variable posX
	VAriable posY
	
	NVar storageHoughResolution
	
	NVar intervalleHoughSlopes
	NVar intervalleHoughIntercepts
	
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	intervalleHoughSlopes = intervalleHoughSlopes/2
	intervalleHoughIntercepts = intervalleHoughIntercepts/2
				
	borneHoughSlopes1 = borneHoughSlopes1 + ((posX - 80)*(intervalleHoughSlopes*2)/storageHoughResolution) - (intervalleHoughSlopes/2)
	borneHoughSlopes2 = borneHoughSlopes1 + ((posX - 80)*(intervalleHoughSlopes*2)/storageHoughResolution) + (intervalleHoughSlopes/2)
				
	borneHoughIntercepts1 = borneHoughIntercepts2 - ((posY - 20)*(intervalleHoughIntercepts*2)/storageHoughResolution) - (intervalleHoughIntercepts/2)
	borneHoughIntercepts2 = borneHoughIntercepts2 - ((posY - 20)*(intervalleHoughIntercepts*2)/storageHoughResolution) + (intervalleHoughIntercepts/2)
				
	miseAjourTransfo()
End // zoomPar2(posX,posY)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// dezoomPar2(posX,posY)
//	Recalcule la transformée de Hough pour avoir un dezoom par deux,
//	zone centrée sur posX posY deux fois plus large
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function dezoomPar2(posX,posY)
	Variable posX
	VAriable posY
	
	NVar storageHoughResolution
	
	NVar intervalleHoughSlopes
	NVar intervalleHoughIntercepts
	
	NVar borneHoughSlopes1
	NVar borneHoughSlopes2
	NVar borneHoughIntercepts1
	NVar borneHoughIntercepts2
	
	intervalleHoughSlopes = intervalleHoughSlopes*2
	intervalleHoughIntercepts = intervalleHoughIntercepts*2
				
	borneHoughSlopes1 = borneHoughSlopes1 + ((posX - 80)*(intervalleHoughSlopes/2)/storageHoughResolution) - (intervalleHoughSlopes/2)
	borneHoughSlopes2 = borneHoughSlopes1 + ((posX - 80)*(intervalleHoughSlopes/2)/storageHoughResolution) + (intervalleHoughSlopes/2)
				
	borneHoughIntercepts1 = borneHoughIntercepts2 - ((posY - 20)*(intervalleHoughIntercepts/2)/storageHoughResolution) - (intervalleHoughIntercepts/2)
	borneHoughIntercepts2 = borneHoughIntercepts2 - ((posY - 20)*(intervalleHoughIntercepts/2)/storageHoughResolution) + (intervalleHoughIntercepts/2)
				
	miseAjourTransfo()
End // dezoomPar2(posX,posY)

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Fonctions de représentations
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// changeStyleGraph(nomBouton, variableChangee)
//	Une des CheckBox controlant les options d'affichage de l'image de
//	la transformée, il faut donc recalculer les options d'affichage en
//	fonction des checkboxes
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function changeStyleGraph(nomBouton, variableChangee) : CheckBoxControl
	String nomBouton
	Variable variableChangee
	
	NVar storageBoolHoughCouleur
	NVar storageBoolHoughReverse
	NVar storageBoolHoughLog
	
	If(storageBoolHoughCouleur==1)
		ModifyImage /W=etudeTransformationHough#imageResultatHough resultatHough2 ctab= {storageBoolHoughLog,*,BlueGreenOrange,storageBoolHoughReverse}, log=storageBoolHoughLog
	Else
		ModifyImage /W=etudeTransformationHough#imageResultatHough resultatHough2 ctab= {storageBoolHoughLog,*,Grays,storageBoolHoughReverse}, log=storageBoolHoughLog
	endIf	
End // changeStyleGraph(nomBouton, variableChangee)

////////////////////////////////////////////////////////////////////////////////////////////////////////
// tracerRosedesPentes(nomBouton)
//	Créer, calcule les waves necessaire à la rose des pentes avant de 
//	la représentée. Elle est composé d'une rose basée sur l'ensemble
//	des données, une seconde rose basée sur la selection en cours et
//	enfin des différents facteurs référence issu du molbag sélectioné
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function tracerRosedesPentes(nomBouton) : ButtonControl
	String nomBouton
	
	NVar storageHoughResolution
	
	Wave resultatHough2Initial
	Wave histoHoughSlopes
	Wave histoHoughIntercepts
	Wave wave0stripped
	
	If(winType("RoseDesPentes") != 0)
		Killwindow RoseDesPentes
	endIf
	
	// Fabrique histo de l'etat initial
	Make /O /N=(DimSize(resultatHough2Initial,0)) histoHoughSlopeRoseGen
	Matrixop/O histoHoughSlopeRoseGen=sumRows(resultatHough2Initial)
	
	Make/O/N=2 SlopeCursorY=0, SlopeCursorX=0
	Make /O /N=(numpnts(histoHoughSlopeRoseGen)) rosePentesX,xplode
	Make /O /N=(numpnts(histoHoughSlopeRoseGen)) rosePentesY,yplode
	
	Duplicate /O histoHoughSlopeRoseGen histoHoughRoseGenSlope, histoHoughRoseGenRayon
	histoHoughRoseGenSlope = -0.01 + (x*(0.02/storageHoughResolution))
	histoHoughRoseGenRayon = sqrt(histoHoughRoseGenSlope*histoHoughRoseGenSlope + histoHoughSlopeRoseGen*histoHoughSlopeRoseGen)
	
	rosePentesX = histoHoughRoseGenRayon * cos(atan(histoHoughRoseGenSlope))
	rosePentesY = histoHoughRoseGenRayon * sin(atan(histoHoughRoseGenSlope))
	
	// Fabrique Cauchy
	CurveFit/Q /W=2 lor, histoHoughSlopeRoseGen/D
	Wave W_coef
	Make/O/N=(numpnts(histoHoughSlopeRoseGen)) cauchyX, cauchyY, cauchyxplode, cauchyyplode
	cauchyx  =histoHoughRoseGenSlope
	Variable k
	For(k=0 ; k<numpnts(cauchyy) ; k+=1)
		cauchyy[k] = W_coef[0] + W_coef[1] / ( (k-W_coef[2])^2 + W_coef[3])
	endFor
	
	cauchyxplode=cauchyy*cos(atan(cauchyx))
	cauchyyplode=cauchyy*sin(atan(cauchyx))
	
	xplode=max(0,(histoHoughRoseGenRayon-cauchyy))*cos(atan(histoHoughRoseGenSlope))
	yplode=max(0,(histoHoughRoseGenRayon-cauchyy))*sin(atan(histoHoughRoseGenSlope))
	Variable maxx=max(wavemax(xplode), abs(wavemin(xplode)))
	Variable maxy=max(wavemax(yplode), abs(wavemin(yplode)))
	Variable contrainte=min((wavemax(wave0stripped)-wavemin(wave0stripped))/maxx,0.5/maxy)
	xplode*=contrainte
	yplode*=contrainte
	cauchyxplode*=contrainte
	cauchyyplode*=contrainte
	rosePentesX*=contrainte
	rosePentesY*=contrainte
	Wave wave0stripped
	//xplode+=wavemin(wave0stripped)
	//cauchyxplode+=wavemin(wave0stripped)
	//rosePentesX+=wavemin(wave0stripped)

	Display /W=(50,50,600,600) /K=1/N=RoseDesPentes rosePentesY vs rosePentesX
	// ajoute le fit
	AppendToGraph /W=RoseDesPentes /C=(0,0,65000) cauchyyplode vs cauchyxplode
	AppendToGraph /W=RoseDesPentes SlopeCursorY vs SlopeCursorX
	ModifyGraph mode(SlopeCursorY)=4,marker(SlopeCursorY)=19;DelayUpdate
	ModifyGraph rgb(SlopeCursorY)=(0,64768,0),useMrkStrokeRGB(SlopeCursorY)=1
	
	// Affiche selection si il y en a une
	GetMarquee /W=etudeTransformationHough#imageResultatHough left,bottom
	//Ajoutte rose pente selection
	Make /O /N=(numpnts(histoHoughSlopes)) rosePentesXSelect
	Make /O /N=(numpnts(histoHoughSlopes)) rosePentesYSelect
	AppendToGraph /W=RoseDesPentes /C=(0,65000,0) rosePentesYSelect vs rosePentesXSelect
	If(V_Flag!=0)
		NVar borneHoughSlopes1
		NVar intervalleHoughSlopes
		NVar storageHoughResolution

		Duplicate /O histoHoughSlopes histoHoughRoseSlopeSelect histoHoughRoseGenRSelect
		Variable i
		For(i=0 ; i<numpnts(histoHoughRoseSlopeSelect) ; i +=1)
			histoHoughRoseSlopeSelect[i] = borneHoughSlopes1 + ((V_left + i) * intervalleHoughSlopes/storageHoughResolution)
			histoHoughRoseGenRSelect[i] = sqrt(histoHoughRoseSlopeSelect[i]*histoHoughRoseSlopeSelect[i] + histoHoughSlopes[i]*histoHoughSlopes[i])
		endFor
	
		rosePentesXSelect = histoHoughRoseGenRSelect * cos(atan(histoHoughRoseSlopeSelect))
		rosePentesYSelect = histoHoughRoseGenRSelect * sin(atan(histoHoughRoseSlopeSelect))
		
		rosePentesXSelect*=contrainte
		rosePentesYSelect*=contrainte
		//rosePentesXSelect+=wavemin(wave0stripped)
	endIf
	
	//Affichage des facteurs d'attribution
	// On ajoutte une courbe par facteur -> ceci est LENT.
	// Mais cela permet d'utiliser un tag qui à la propriété de suivre une trace.
	// Ainsi on a un label qui est collé à la courbe dont il donne le nom, 
	// cela rend possible le placement des labels (qui sont des tag).
	// Le défaut majeur est que ce systeme est lent et que lors du chargement d'un
	//nouveau molbag de facteur c'est franchement ridicule que quelquechose d'aussi
	// simple soit aussi lent.
	Wave listeMasseFacteursRosePentes
	Wave/T listeFacteursRosePentes
	Variable l = 0
	For(l=0 ; l<numpnts(listeMasseFacteursRosePentes) ; l+=1)
		String nomWaveX = "roseDesPentesFacteur" + num2str(l) + "X"
		String nomWaveY = "roseDesPentesFacteur" + num2str(l) + "Y"
		Make /O /N=2 $nomWaveX
		Make /O /N=2 $nomWaveY
		Wave roseEnConstructionX  = $nomWaveX
		Wave roseEnConstructionY  = $nomWaveY
		Variable penteEnconstruction = (listeMasseFacteursRosePentes[l] - round(listeMasseFacteursRosePentes[l])) / listeMasseFacteursRosePentes[l]
		FindLevel /P /Q cauchyx, penteEnconstruction
		roseEnConstructionX[0]=0//+wavemin(wave0stripped)
		roseEnConstructionX[1]=rosepentesX[V_LevelX]-wavemin(wave0stripped)
		roseEnConstructionY[0]=0
		roseEnConstructionY[1]=rosepentesY[V_LevelX]
		
		AppendToGraph /W=RoseDesPentes /C=(0,0,0) roseEnConstructionY vs roseEnConstructionX
		
		Tag /W=RoseDesPentes /N=$nomWaveY $nomWaveY,cauchyxplode[V_LevelX] ,listeFacteursRosePentes[l]+"\M,"+num2str(penteEnconstruction)
	endFor
End //tracerRosedesPentes(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fonctions Autres
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// boolSelectionChanged()
//	Fonction permetant de suivre si la selection a été changée. Pour ce
//	faire elle compare à la valeur en memoire dans les variables globales
//	memoireSelectHough. Si la selection a changée alors la fonction
//	place la nouvelle position en memoire et retourne 1, sinon elle ne
//	change rien et renvoit 0.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function boolSelectionChanged()
	NVar memoireSelectHoughTop
	NVar memoireSelectHoughBottom
	NVar memoireSelectHoughLeft
	NVar memoireSelectHoughRight
	
	GetMarquee /W=etudeTransformationHough#imageResultatHough left,bottom
	If(V_Flag==0)
		// Pas de selection
		memoireSelectHoughTop = -1
		memoireSelectHoughBottom = -1
		memoireSelectHoughLeft = -1
		memoireSelectHoughRight = -1
	Else
		If(memoireSelectHoughTop!=V_top || memoireSelectHoughBottom!=V_bottom || memoireSelectHoughLeft!=V_left || memoireSelectHoughRight!=V_right)
			memoireSelectHoughTop = V_top
			memoireSelectHoughBottom = V_bottom
			memoireSelectHoughLeft = V_left
			memoireSelectHoughRight = V_right
			
			return 1
		Else
			return 0
		endIf
	endIf
End // boolSelectionChanged()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// chargerMolbagFacteurAff(nomMenu,indexSelect,strSelect)
//	Fonction qui prépare les fonctions d'affichage des slopes référence
//	dans la rose des pentes à partir d'un molbag. Pour ce faire elle va
//	calculer la masse des différentes molécules et en déduire leur pente.
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
// nettoyageHough()
//	Fonctions de nettoyage des différentes waves utilisées
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function nettoyageHough()
	// NETTOYAGE
	Variable indexCompteur = 0
	
	Variable boolTestPresenceFonction = 1
	Do
		String nomWaveASupprimer1 =  "roseDesPentesFacteur" + num2str(indexCompteur) + "X"
		String nomWaveASupprimer2 =  "roseDesPentesFacteur" + num2str(indexCompteur) + "Y"
		boolTestPresenceFonction = waveexists($nomWaveASupprimer2)
		
		If(boolTestPresenceFonction !=0)		
			KillWaves /Z $nomWaveASupprimer1 $nomWaveASupprimer2		
			indexCompteur+=1
		endIf
	While(boolTestPresenceFonction != 0)

	KillWindow etudeTransformationHough#imageResultatHough
	KillWindow etudeTransformationHough#histoInterceps
	KillWindow etudeTransformationHough#histoSlopes

	KillWaves /Z resultatHough2,rosePentesYSelect,rosePentesXSelect, rosePentesX,xplode,conteneurTicksLabelsIntercepts,conteneurTicksSlopes,listeFacteursRosePentes,listeMasseFacteursRosePentes
	KillWaves /Z rosePentesY,yplode, cauchyX, cauchyY, cauchyxplode, cauchyyplode,histoHoughSlopeRoseGen,memoireTestHough,histoHoughAffichageHoriz2,intercepts3,resultatHough2Initial,SlopeCursorY, SlopeCursorX
	KillWaves /Z histoHoughAffichageHoriz1,histoHoughAffichageVert2,histoHoughAffichageVert1,histoHoughSlopes,histoHoughIntercepts,conteneurTicksIntercepts,histoHoughRoseGenRayon,histoHoughRoseGenSlope,resultatHough2Selection,slopes3
	KillWaves /Z conteneurTicksLabelsSlopes,conteneurTicksInterceptsHisto,conteneurTicksLabelsInterHisto,conteneurTicksSlopesHisto,conteneurTicksLabelsSlopesHisto,histoHoughRoseSlopeSelect,histoHoughRoseGenRSelect
	KillWaves /Z intermarks, slomarks, listeSlopesFacteursHT 
End // nettoyageHough()

Function HoughTcaller(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			etudeTransfoHough()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Window etudeTransformationHough() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(29,88,1279,798)
	ShowTools/A
	GroupBox conteneurParaHough,pos={710,5},size={130,100},title="Hough parameters"
	GroupBox conteneurParaHough,fStyle=1
	TitleBox titreResolutioHough,pos={715,25},size={80,13},title="Hough resolution"
	TitleBox titreResolutioHough,frame=0
	SetVariable setResolution,pos={715,45},size={120,16},title=" "
	SetVariable setResolution,limits={50,10000,50},value= storageHoughResolution,proc=changeHTreso
	TitleBox titreResolutioMassHough,pos={715,65},size={73,13},title="Mass resolution"
	TitleBox titreResolutioMassHough,frame=0
	ValDisplay affResoEnMasse,pos={708,85},size={127,14},bodyWidth=120,title=" "
	ValDisplay affResoEnMasse,limits={0,0,0},barmisc={0,1000},value= #"resoEnMasse"
	GroupBox conteneurZOOM,pos={710,115},size={130,100},title="Tools ZOOM",fStyle=1
	Button boutonZOOM,pos={715,135},size={120,25},proc=zoomSelection,title="ZOOM Selection"
	Button boutonZOOM,fColor=(20000,45000,3000)
	Button boutonDEZOOM,pos={715,160},size={120,25},proc=dezoomSelection,title="deZOOM Selection"
	Button boutonDEZOOM,fColor=(65000,45000,3000)
	Button boutonREINIT,pos={715,185},size={120,25},proc=zoomReinit,title="ZOOM ReInitialisation"
	Button boutonREINIT,fColor=(40000,20000,30000)
	GroupBox conteneurVisu,pos={710,225},size={130,140},title="Vizualisation"
	GroupBox conteneurVisu,fStyle=1
	CheckBox checkCouleur,pos={715,247},size={42,14},proc=changeStyleGraph,title="Color"
	CheckBox checkCouleur,variable= storageBoolHoughCouleur
	CheckBox checkReverse,pos={715,265},size={88,14},proc=changeStyleGraph,title="Inverted colors"
	CheckBox checkReverse,variable= storageBoolHoughReverse
	CheckBox checkLog,pos={715,283},size={66,14},proc=changeStyleGraph,title="Log Scale"
	CheckBox checkLog,variable= storageBoolHoughLog
	Button boutonRosePentes,pos={715,300},size={120,30},proc=tracerRosedesPentes,title="Rose des Pentes"
	Button boutonRosePentes,fColor=(20000,45000,23000)
	PopupMenu choixMolBag,pos={709,335},size={126,21},bodyWidth=120,proc=chargerMolbagFacteurAff,title=" "
	PopupMenu choixMolBag,mode=1,popvalue="Facteurs en Memoire",value= #"\"Facteurs en Memoire;Molbag_Current;Molbag_ClassicOrga;Molbag_Factors;Molbag_CHNO;Molbag_ManualAttrib1;Molbag_CalibCOSI;Molbag_ManualCH;Molbag_Last_Attribution;Molbag_ManualCHNa;Molbag_Attrib2;Molbag_attribCH;Molbag_attribCHN;Molbag_certitude103;Molbag_MgAlCaFe;\""
	GroupBox conteneurPos,pos={845,5},size={400,40},title="Position curseur"
	GroupBox conteneurPos,fStyle=1
	ValDisplay affPosYIntercept,pos={881,25},size={149,14},bodyWidth=90,title="Intercept "
	ValDisplay affPosYIntercept,fStyle=1,limits={0,0,0},barmisc={0,1000}
	ValDisplay affPosYIntercept,value= #"posYIntercept"
	ValDisplay affPosXSlope,pos={1090,25},size={130,14},bodyWidth=90,title="Slope "
	ValDisplay affPosXSlope,fStyle=1,limits={0,0,0},barmisc={0,1000}
	ValDisplay affPosXSlope,value= #"posXSlope"
	SetWindow kwTopWin,hook(imageResultatHough)=hookImageResultatHough
	Display/W=(5,5,705,705)/HOST=#  histoHoughAffichageVert2 vs histoHoughAffichageVert1
	AppendToGraph histoHoughAffichageHoriz1 vs histoHoughAffichageHoriz2
	AppendToGraph plotinterdmvm vs plotslodmvm
	AppendImage resultatHough2
	ModifyImage resultatHough2 ctab= {0,*,Grays,1}
	ModifyGraph userticks(left)={conteneurTicksIntercepts,conteneurTicksLabelsIntercepts}
	ModifyGraph userticks(bottom)={conteneurTicksSlopes,conteneurTicksLabelsSlopes}
	ModifyGraph margin(left)=80,margin(bottom)=80,margin(top)=20,margin(right)=20
	ModifyGraph mode(plotinterdmvm)=3
	ModifyGraph marker(plotinterdmvm)=39
	ModifyGraph rgb(histoHoughAffichageVert2)=(15000,50000,15000),rgb(histoHoughAffichageHoriz1)=(15000,50000,15000)
	ModifyGraph mirror=2
	ModifyGraph fSize=10
	Label left "\\Z14Intercepts (SI)"
	Label bottom "\\Z14Slopes (SI)"
	RenameWindow #,imageResultatHough
	SetActiveSubwindow ##
	Display/W=(845,49,1245,377)/HOST=#  histoHoughIntercepts
	ModifyGraph userticks(bottom)={conteneurTicksInterceptsHisto,conteneurTicksLabelsInterHisto}
	ModifyGraph margin(left)=30,margin(right)=5
	ModifyGraph mode=5
	ModifyGraph notation(left)=1
	ModifyGraph btLen(left)=2
	TextBox/C/N=text0/F=0/B=1/A=MT "Histo Intercepts (vertical)"
	RenameWindow #,histoInterceps
	SetActiveSubwindow ##
	Display/W=(845,377,1245,705)/HOST=#  histoHoughSlopes
	ModifyGraph userticks(bottom)={conteneurTicksSlopesHisto,conteneurTicksLabelsSlopesHisto}
	ModifyGraph margin(left)=30,margin(right)=5
	ModifyGraph mode=5
	ModifyGraph notation(left)=1
	ModifyGraph btLen(left)=2
	TextBox/C/N=text0/F=0/B=1/A=MT "Histo Slopes (horizontal)"
	RenameWindow #,histoSlopes
	SetActiveSubwindow ##
EndMacro

Function PopMenuProc_2(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			molbag2perm(popstr)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//////////
//rework du 14 fev 2019
/////
function MolList4Attributeur(quant,bag)
variable quant
string bag
string sciblebag="Ciblebag_"+bag
wave ciblebag=$sciblebag
InitMolList(quant,bag)
string lasel="selmol"+num2str(quant)
string lalist="mollist"+num2str(quant)
string latitle="moltitles"+num2str(quant)
wave selmol=$lasel
wave/T mollist=$lalist
wave/T moltitles=$latitle
insertpoints numpnts(moltitles),2, moltitles
moltitles[numpnts(moltitles)-2]="Minimum"
moltitles[numpnts(moltitles)-1]="Maximum"
insertpoints/M=1 dimsize(selmol,1),2, selmol,mollist
selmol[][dimsize(selmol,1)-2]=2
selmol[][dimsize(selmol,1)-1]=2
Duplicate/O/FREE bag2perm("bag"), perm
selmol[perm[dimsize(selmol,0)-1][1]][dimsize(selmol,1)-1]=0
if(dimsize(ciblebag,1)==5)
	insertpoints/M=1 5,2,ciblebag
	ciblebag[][5]=0
	ciblebag[][6]=inf
endif
if(dimsize(ciblebag,1)==7)
	mollist[][dimsize(selmol,1)-2]=num2str(ciblebag[x][5])
	mollist[][dimsize(selmol,1)-1]=num2str(ciblebag[x][6])
endif
end
///////

function getHoughMatch()
	
	NVar posXslope
	NVar posYintercept
	NVar precinter
	NVar precslo
	NVar inco

variable k,n,t=ticks
wave slopes, intercepts, modules,centremasse, wave0stripped, wave1stripped, dataindex
duplicate/O slopes match
match= abs(slopes-posXslope)<precslo && abs(intercepts-posYintercept-inco*(posXslope-slopes))<precinter
Matrixop/O checkslo=(sumcols(match))^t
//simplifysmatch()
Make/D/O/N=(numpnts(modules)) modtrans, pentrans, interceptrans, centremassetrans, xxtrans
pentrans=match[x]*slopes[x]
modtrans=match[x]*modules[x]
interceptrans=match[x]*intercepts[x]
centremassetrans=match[x]*centremasse[x]
xxtrans=match[x]*x
sort modtrans, modtrans, pentrans, interceptrans, xxtrans, centremassetrans
deletepoints 0,numpnts(modules)-sum(checkslo), modtrans, pentrans, interceptrans, xxtrans, centremassetrans
//sort/R modtrans, modtrans, pentrans, interceptrans
Make/O/N=0 roi0,roi1,roi2,roipnts
n=numpnts(checkslo)
	for(k=0;k<n;k+=1)
		if(checkslo[k]!=0)
		insertpoints 0,1,roi0,roi1,roi2,roipnts
		roipnts[0]=dataindex[k]
		roi0[0]=wave0stripped[dataindex[k]]
		roi1[0]=wave1stripped[dataindex[k]]
		endif
	endfor
	sort roi0, roi0,roi1
	roi2=roi0-round(roi0)
//killwaves match
end

function listelesdonnees()
//list the spectrum data
string lesuns=wavelist("*_1stripped",";",""),laun
variable k, n=itemsinlist(lesuns), taille
Make/O/T/N=(n) liste_donnees
wave/T liste_donnees
make/O/T/N=6 seladvtitles={"","Name","Points","Spe","DM",""}
make/T/O/N=(n,6) advdatalist=""
make/O/N=(n,6) seladvdatalist=0
	for(k=0;k<n;k+=1)
		laun=stringfromlist(k,lesuns)
		advdatalist[k][2]=num2str(numpnts($laun))
		taille=strlen(laun)-11
		laun=laun[0,taille]
		liste_donnees[k]=laun
		advdatalist[k][1]=laun
	endfor
seladvdatalist[][0]=64
seladvdatalist[][3]=32+16*(strsearch(TraceNameList("agregateur", ";", 1), advdatalist[x][1], 0)>0)
seladvdatalist[][4]=32+16*(strsearch(TraceNameList("dmvm", ";", 1), advdatalist[x][1], 0)>0)
seladvdatalist[][5]=32
//list the chromatogram data
string leschro=wavelist("*_1chroma",";",""),lechro
string name,stiti,smama,spopo,spipi
n=itemsinlist(leschro)
Make/O/T/N=(n) liste_chroma
wave/T liste_chroma
make/O/T/N=8 chromatitles={"","Name","Time","Mass","Points","Peaks","",""}
make/T/O/N=(n,8) advchromalist=""
make/O/N=(n,8) advchromasel=0
	for(k=0;k<n;k+=1)
		name=ReplaceString("_1chroma",stringfromlist(k,leschro),"")
		advchromalist[k][1]=name
		stiti=name+"_3chroma"
		wave titi =$stiti
		smama=name+"_4chroma"
		wave mama=$smama
		spopo=name+"_6chroma"
		wave popo=$spopo
		spipi=name+"_11chroma"
		wave pipi=$spipi
		if(waveexists(titi))
			advchromalist[k][2]=num2str(wavemin(titi))+"-"+num2str(wavemax(titi))
		else
			advchromalist[k][2]="N/A"
		endif
		if(waveexists(mama))
			advchromalist[k][3]=num2str(wavemin(mama))+"-"+num2str(wavemax(mama))
		else
			advchromalist[k][3]="N/A"
		endif
		if(waveexists(popo))
			advchromalist[k][4]=num2str(dimsize(popo,0))
		else
			advchromalist[k][4]="N/A"
		endif
		if(waveexists(popo))
			advchromalist[k][5]=num2str(dimsize(pipi,0))
		else
			advchromalist[k][5]="N/A"
		endif
		endfor
advchromasel[][0]=64
end

Function InitAdvDataMan()
newpanel/K=2/N=AdvancedManager
Movewindow/W=AdvancedManager 1083,  567.5,  1389,  704
getwindow AdvancedManager, wsizeDC
variable width=V_right-V_left, height=V_bottom-V_top
TabControl tab0 win=AdvancedManager, pos={0,0}, size={width,height}, tabLabel(0)="Spectra", tabLabel(1)="Chromatograms", tabLabel(2)="ROIs",proc=TabAdvDataMan
setwindow AdvancedManager, hook(majdataman)=MajAdvDataMan
SetControlsDATA()
end

Function RazcontrolsDATA()
string lescontrols=ControlNameList("AdvancedManager",";", "TabCon_*DATA")
variable k, n=itemsinlist(lescontrols)
string lecontrol
for(k=0;k<n;k+=1)
    lecontrol=stringfromlist(k,lescontrols)
    killcontrol/W=AdvancedManager $lecontrol
endfor
end

Function SetControlsDATA()
controlinfo/W=AdvancedManager tab0
variable width=V_width, height=V_height, tab=V_value,nbbutt
RazcontrolsDATA()
    switch(tab)
        case 0:
        Button TabCon_Load_DATA win=AdvancedManager, pos={2,20}, size={60,(height-23)/5}, title="Load",proc=ButtonLoader
        Button TabCon_Activate_DATA win=AdvancedManager, pos={2,20+(height-22)/5}, size={60,(height-23)/5}, title="Activate",proc=ActiveLaDataButton,fColor=(65280,43520,0)
        Button TabCon_ToROI_DATA win=AdvancedManager, pos={2,20+2*(height-23)/5}, size={60,(height-23)/5}, title="To ROI",proc=loadInRoi
        Button TabCon_Duplicate_DATA win=AdvancedManager, pos={2,20+3*(height-23)/5}, size={60,(height-23)/5}, title="Save active",proc=ButtonDupeRname
        Button TabCon_Kill_DATA win=AdvancedManager, pos={2,20+4*(height-23)/5}, size={60,(height-23)/5}, title="Kill",proc=ButtonKill,fColor=(65280,0,0)
        ListBox TabCon_List_DATA win=AdvancedManager, pos={62,20}, size={width-63,(height-22)}, mode=1, listwave=advdatalist,selwave=seladvdatalist,titlewave=seladvtitles,widths={22,140,71,22,22,22},proc=ListBoxProcDataMan
            break
        case 1:
        nbbutt=5
        Button TabCon_Load_CHROMA win=AdvancedManager, pos={2,20}, size={60,(height-23)/nbbutt}, title="Load",proc=ButtonChromaLoader
        Button TabCon_Activate_CHROMA win=AdvancedManager, pos={2,20+(height-22)/nbbutt}, size={60,(height-23)/nbbutt}, title="Build Map",proc=ButtonChromaMap
        Button TabCon_ToROI_CHROMA win=AdvancedManager, pos={2,20+2*(height-23)/nbbutt}, size={60,(height-23)/nbbutt}, title="Timeless\rRoi",proc=ButtonShowFittedPeaks
        Button TabCon_Duplicate_CHROMA win=AdvancedManager, pos={2,20+3*(height-23)/nbbutt}, size={60,(height-23)/nbbutt}, title="TBD"//,proc=
        Button TabCon_Kill_CHROMA win=AdvancedManager, pos={2,20+4*(height-23)/nbbutt}, size={60,(height-23)/nbbutt}, title="Kill",proc=ButtonChromaKill,fColor=(65280,0,0)
        ListBox TabCon_List_CHROMA win=AdvancedManager, pos={62,20}, size={width-63,(height-22)}, mode=1, listwave=advchromalist,selwave=advchromasel,titlewave=chromatitles,widths={22,140,71,71,71,22,22,22}//,proc=
            break
        case 2:
        //Button
            break
    endswitch
end

Function ButtonChromaMap(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
		controlinfo/W=AdvancedManager TabCon_List_CHROMA
		wave/t lasel=$S_value
		string name=lasel[V_value][1]
		InitChromaMapBuilder(name)
			break
	endswitch
	return 0
End

Function TabAdvDataMan(tca) : TabControl
    STRUCT WMTabControlAction &tca
    switch( tca.eventCode )
        case 2: // mouse up
            Variable tab = tca.tab
            SetControlsDATA()
            break
    endswitch
    return 0
End

Function MajAdvDataMan(s)
    STRUCT WMWinHookStruct &s
    Variable hookResult = 0
    switch(s.eventCode)
        case 0:
            break
        case 6:
            getwindow AdvancedManager, wsizeDC
            variable width=V_right-V_left, height=V_bottom-V_top
            TabControl tab0 win=$s.winname, pos={0,0}, size={width,height}, tabLabel(0)="Spectra", tabLabel(1)="Chromatograms", tabLabel(2)="ROIs",proc=TabAdvDataMan
            SetControlsDATA()
            break
    endswitch
    return 0
End

function resortDATAbycol(selwave,listwave,col)
variable col
wave/T listwave
wave selwave
variable n=dimsize(listwave,0),m
Make/FREE/O/D/N=(n) index
if(numtype(str2num(listwave[0][col]))==0)
	Make/FREE/O/D/N=(n) key
	key=str2num(listwave[x][col])
	Makeindex key, index
	m=numpnts(key)
else
	Make/FREE/O/T/N=(n) keyt
	keyt=listwave[x][col]+num2str(selwave[x][col])
	Makeindex/A keyt, index
	m=numpnts(keyt)
endif
Duplicate/FREE/O/T listwave newlist
Duplicate/FREE/O selwave newsel
if(index[0]<index[m-1])
    wavetransform/O flip index
endif
newsel=selwave[index[x]][y]
newlist=listwave[index[x]][y]
duplicate/O/T newlist listwave
Duplicate/O newsel selwave
end

Window advancedmanager() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1443,757,1852,939)
	ListBox list0,pos={1,1},size={324,181},proc=ListBoxProcDataMan
	ListBox list0,listWave=root:advdatalist,selWave=root:seladvdatalist
	ListBox list0,titleWave=root:seladvtitles,row= 130,mode= 1,selRow= 132
	ListBox list0,widths={22,140,71,24,25,24},userColumnResize= 1
	Button buttonLoadData,pos={328,1},size={80,25},proc=ButtonLoader,title="Load"
	Button buttonLoadData,help={"Load an XY text file with mass and intensities."}
	Button buttonLoadData,fColor=(32768,65280,0)
	Button buttonFatNoise1,pos={328,28},size={80,25},proc=ActiveLaDataButton,title="Activate"
	Button buttonFatNoise1,help={"Remove the less intense points so the new intensity distribution is strictly convex"}
	Button buttonFatNoise1,fColor=(65280,43520,0)
	Button buttonLoadinROI,pos={328,55},size={80,25},proc=loadInRoi,title="To ROI"
	Button buttonLoadinROI,help={"Load an XY text file with mass and intensities."}
	Button buttonLoadinROI,fColor=(32768,65280,0)
	Button buttonRename,pos={328,82},size={80,25},proc=ButtonDupeRname,title="Dupe active"
	Button buttonRename,help={"Generates an other sample data from the current data."}
	Button buttonRename,fColor=(32768,65280,0)
	Button buttonKiller,pos={328,109},size={80,25},proc=ButtonKill,title="Kill"
	Button buttonKiller,help={"Removes the selected sample data."}
	Button buttonKiller,fColor=(65280,0,0)
EndMacro



Function switchpeakstick(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
controlinfo/W=attributeur checkFindInData
nvar findindata
variable/G spectrummode
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			switch(popnum)
				case 1:
					ModifyGraph/W=agregateur mode(wave1stripped)=1
					spectrummode=1
					break
				case 2:
					ModifyGraph/W=agregateur mode(wave1stripped)=0
					spectrummode=0
					findindata*=2
					break
			endswitch
			break
		case -1:
			break
	endswitch

	return 0
End

Function fish(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = exp(-l)*l^x/factorial(round(x))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 1
	//CurveFitDialog/ w[0] = l

	return exp(-w[0])*w[0]^x/factorial(round(x))
End


Function ListBoxProcAgregToElab(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
wave/Z Molecule,Index_formules,incox,incoy,incoxdm, wave1stripped, wave0stripped
string lawave
nvar inco

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listwave = lba.listWave
	WAVE/Z selwave = lba.selWave

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 3: // double click
			break
		case 4: // cell selection
			controlinfo/W=molmanager list0
			if(v_value>0)
			lawave=stringfromlist(v_value-1,wavelist("Molecule_*",";",""))
			molecule=0
			tuelesspec()
			duplicate/O $lawave molecule
			genestringmol()
			inco=wave0stripped(prochepic(index_formules[v_value-1],0))
			incox=inco
			incoxdm=defmass(incox)
			incoy=wave1stripped(prochepic(inco,0))			
			endif
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 13:
			if(row==0 && col==0)
				selwave[][0]=selwave[0][0]
			endif
			Majlegende()
			break
	endswitch

	return 0
End

Function ButtonNorm(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z wave1stripped
			variable xeuma=wavemax(wave1stripped)
			duplicate/O wave1stripped roi1
			duplicate/O wave0stripped roi0
			roi1/=xeuma
			break
	endswitch

	return 0
End


Window Ticked_cropper() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1168,651,1511,685)
	CheckBox check2,pos={1,1},size={81,14},title="Crop in mass:",value= 1
	CheckBox check3,pos={1,17},size={95,14},title="Crop in intensity:",value= 0
	SetVariable setvar0,pos={84,1},size={96,16},title="from",value= _NUM:0
	SetVariable setvar1,pos={181,1},size={108,16},title="to",value= _NUM:500
	SetVariable setvar2,pos={97,17},size={96,16},title="from",value= _NUM:0
	SetVariable setvar3,pos={194,17},size={96,16},title="to",value= _NUM:0
	Button button0,pos={292,2},size={50,30},proc=ButtonProc_14,title="GO !"
	Button button0,fColor=(65280,0,0)
EndMacro

function tickedcropper()
variable bm1,bm2, bi1, bi2, cropmass, cropint
wave/T advdatalist
wave seladvdatalist
variable n=dimsize(advdatalist,0),k
string lenomx,lenomy
controlinfo/W=ticked_cropper setvar0
bm1=v_value
controlinfo/W=ticked_cropper setvar1
bm2=v_value
controlinfo/W=ticked_cropper setvar2
bi1=v_value
controlinfo/W=ticked_cropper setvar3
bi2=v_value
controlinfo/W=ticked_cropper check2
cropmass=v_value
controlinfo/W=ticked_cropper check3
cropint=v_value
for(k=0;k<n;k+=1)
	if(seladvdatalist[k][5]==48)
		lenomx=advdatalist[k][1]+"_0stripped"
		lenomy=advdatalist[k][1]+"_1stripped"
		if(cropmass)
			cropbetween(bm1,bm2,$lenomy,$lenomx)
			duplicate/O outx $lenomy
			duplicate/O outy $lenomx
		endif
		if(cropint)
			cropbetween(bi1,bi2,$lenomx,$lenomy)
			duplicate/O outx $lenomx
			duplicate/O outy $lenomy
		endif
		killwaves outx, outy
	endif
endfor
end

Function ButtonProc_14(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			tickedcropper()
			listelesdonnees()
			killwindow ticked_cropper
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_15(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			execute "Ticked_cropper()"
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function COSIMAsumTicked()
wave seladvdatalist, transint
wave/T advdatalist
variable n=dimsize(seladvdatalist,0),k, counter=0
string lay, lax
	for(k=0;k<n;k+=1)
		if(seladvdatalist[k][5]==48)
			lay=advdatalist[k][1]+"_1stripped"
			if(counter==0)
				lax=advdatalist[k][1]+"_0stripped"
				Duplicate/O $lax COSIMAsum_0stripped
				make/O/D/N=(numpnts(COSIMAsum_0stripped)) COSIMAsum_1stripped
			endif
			duplicate/O $lay transint
			COSIMAsum_1stripped+=transint
			counter+=1
		endif
	endfor
killwaves transint
duplicate/O COSIMAsum_1stripped COSIMAverage_1stripped
COSIMAverage_1stripped/=counter
duplicate/O COSIMAsum_0stripped COSIMAverage_0stripped
listelesdonnees()
end

function IsoRoiCarb()
nvar inco
wave wave1stripped, wave0stripped,wave2stripped
make/O/D/N=4 roipnts,roi0, roi1, roi2, errorCarb
roipnts={inco, inco+1.003355, inco+12, inco +12+1.003355}
roipnts=prochepic(roipnts,0)
roi1=wave1stripped[roipnts[x]]
roi0=wave0stripped[roipnts[x]]
roi2=wave2stripped[roipnts[x]]
errorcarb={inco, inco+1.003355, inco+12, inco +12+1.003355}
errorcarb=abs(1e6*(roi0-errorcarb)/errorcarb)
end

function CarbRoi2Ratio()
wave roi0, roi1, roi2
variable ratio1312
ratio1312=roi1[3]/roi1[2] - roi1[1]/roi1[0]
return ratio1312
end

function CarbRoi2MeanInt()
wave roi0, roi1, roi2
variable meanInt
meanint=mean(roi1)
return meanint
end

function CarbRoi2MeanErr()
wave errorcarb
variable meanErr
meanErr=mean(errorcarb)
return meanErr
end

function CarbRoi2MeanMass()
wave roi0
variable meanMass
meanmass=mean(roi0)
return meanmass
end

function spanCarbRatio()
wave dataindex, wave0stripped
make/O/D/N=0 ratio1312, meanint, meanErr, meanMass
//display ratio1312 vs meanint
//display ratio1312 vs meanErr
//newmovie/A/F=10
variable n=numpnts(dataindex),k
nvar inco
for(k=0;k<n;k+=1)
	inco=wave0stripped[dataindex[k]]
	insertpoints 0,1, ratio1312, meanint, meanErr, meanmass
	IsoRoiCarb()
	ratio1312[0]=CarbRoi2Ratio()
	meanint[0]=CarbRoi2MeanInt()
	meanErr[0]=CarbRoi2MeanErr()
	meanmass[0]=CarbRoi2MeanMass()
	//doupdate
	//addmovieframe
endfor
end

function cropaccuratio(seuil)
variable seuil
wave ratio1312, meanint, meanErr, meanmass
sort meanErr, ratio1312, meanint, meanErr, meanmass
findlevel/Q meanerr, seuil
deletepoints v_levelx, numpnts(ratio1312), ratio1312, meanint, meanErr, meanmass
end

function funnel(chance,fac)
variable chance,fac
wave ratio1312,meanint, meanErr, meanmass
sort meanint,ratio1312,meanint, meanErr, meanmass
duplicate/O ratio1312, funnelp,funneln, testfunnel
funnelp=chance+2*sqrt((1-chance)*chance)/sqrt(meanint/fac)
funneln=chance-2*sqrt((1-chance)*chance)/sqrt(meanint/fac)
testfunnel=(ratio1312<funnelp) && (ratio1312>funneln)
return sum(testfunnel)/numpnts(testfunnel)
end

function dirstreak2ROI()
nvar inco
wave currentseg, wave0stripped, wave1stripped,wave2stripped, roi1, roi0, roi2, treeerror, treeindata
variable deltam=currentseg[0][1]-currentseg[1][1], mcounter=0,maxi=wavemax(wave0stripped), mini=wavemin(wave0stripped), target
variable referror
findvalue/V=(inco+deltam) treeindata
referror=treeerror(v_value)
for(mcounter=inco; mcounter<maxi && mcounter>mini; mcounter+=deltam)
	target=wave0stripped[prochepic(mcounter,pi)]
	if(abs(abs(1e6*(mcounter-target)/mcounter)-referror)/referror <5) ///ARBITRARY PARAMETER
	insertpoints 0,1, roi1, roi0, roi2
	roi0[0]=target
	roi1[0]=wave1stripped[prochepic(mcounter,pi)]
	roi2[0]=wave2stripped[prochepic(mcounter,pi)]
	doupdate
	endif
endfor
for(mcounter=inco-deltam; mcounter<maxi && mcounter>mini; mcounter-=deltam)
	target=wave0stripped[prochepic(mcounter,pi)]
	if(abs(abs(1e6*(mcounter-target)/mcounter)-referror)/referror <5) ///ARBITRARY PARAMETER
	insertpoints 0,1, roi1, roi0, roi2
	roi0[0]=target
	roi1[0]=wave1stripped[prochepic(mcounter,pi)]
	roi2[0]=wave2stripped[prochepic(mcounter,pi)]
	doupdate
	endif
endfor
sort roi0, roi0, roi1, roi2
calcinfo()
end
Function flory(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = b*a^2*(x/d-c)*(1-a)^(x/d-c-1)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = a
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = c
	//CurveFitDialog/ w[3] = d

	return w[1]*w[0]^2*(x/w[3]-w[2])*(1-w[0])^(x/w[3]-w[2]-1)
End

function do2ticked()
wave/T advdatalist
wave seladvdatalist
wave W_coef
variable n=dimsize(advdatalist,0),k
string lenomx,lenomy
for(k=0;k<n;k+=1)
	if(seladvdatalist[k][5]==48)
		lenomx=advdatalist[k][1]+"_0stripped"
		lenomy=advdatalist[k][1]+"_1stripped"
	endif
endfor
end

function/S listticked()
wave/T advdatalist
wave seladvdatalist
wave W_coef
variable n=dimsize(advdatalist,0),k
string lenomx,lenomy, list=""
for(k=0;k<n;k+=1)
	if(seladvdatalist[k][5]==48)
		lenomx=advdatalist[k][1]+"_0stripped"
		lenomy=advdatalist[k][1]+"_1stripped"
		list+=lenomx+";"+lenomy+";"
	endif
endfor
return list
end

function saveticked()
string list=listticked()
save/T/B/W list as "newset.dat"
end

function saveproject()
string list=""
list+=wavelist("*_1stripped",";","")+wavelist("*_0stripped",";","")+wavelist("*bag_*",";","")+wavelist("*chroma",";","")+wavelist("WinLayout",";","")
save/T/B/W list as "newproject.dat"
end

function loadproject()
Loadwave/T/O
listelesdonnees()
end

function histomolposs()
variable bmin=150, bmax=700, pas=1
wave massetot, lesmassescriblees
nvar inco
make/O/D/N=(bmax-bmin) nbposs, massposs
nbposs=Nan
massposs=bmin+x*pas
variable k
	for(k=0;k<((bmax-bmin)/pas);k+=1)	// Initialize variables;continue test
		inco=massposs[k]
		analysemass("arg")
		duplicate/O lesmassescriblees testdesmassescriblees
		testdesmassescriblees= lesmassescriblees<(massposs[k]+1) && lesmassescriblees>massposs[k]
		nbposs[k]=sum(testdesmassescriblees)
		doupdate
	endfor												// Execute body code until continue test is FALSE
end

Function ChangeResoHoughT(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			miseAjourTransfo()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Maxwell_Boltzmann(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = h*x^2*exp(-((x-c)/l)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = l
	//CurveFitDialog/ w[1] = h
	//CurveFitDialog/ w[2] = c

	return w[1]*x^2*exp(-((x-w[2])/w[0])^2)
End

Function PowerFlory(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = A*p^((x-c)/m)*((x-c)/m)^k
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = p
	//CurveFitDialog/ w[2] = c
	//CurveFitDialog/ w[3] = k
	//CurveFitDialog/ w[4] = m

	return w[0]*w[1]^((x-w[2])/w[4])*((x-w[2])/w[4])^w[3]
End


function tree2roi()
variable n=numpnts(tree),k
wave treeindata, treeintens, defmasstree, wave0stripped, wave1stripped
duplicate/O treeindata roi0
duplicate/O treeintens roi1
duplicate/O defmasstree roi2
Make/O/D/N=0 dejavu
for(k=0;k<n;k+=1)
	findvalue/V=(roi0[n-k-1]) dejavu
	if(V_value==-1 && numtype(roi0[n-k-1])==0)
		insertpoints 0,1, dejavu
		dejavu[0]=roi0[n-k-1]
	else
		deletepoints n-k-1,1, roi0,roi1,roi2
	endif
endfor
killwaves dejavu
roi2=defmass(roi0)
sort roi0, roi0, roi1, roi2
end

////////////////////////////////////////////////////////////////////////////////////////////////////////
// fullautoROI
//	Copie de fullauto pour attribuer le contenu de la ROI
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function fullautoROI(roim,roii,roidm,seuil)// ATTENTION POSSIBLES PROBLEMES DU AU FAIT QUE NBPROP PEUT ETRE SUPERIEUR A J, CONTOURNEMENT PAR STR2NUM !!!
	wave roim,roii,roidm
	Variable seuil
	Variable k, n=numpnts(roim), j

	Nvar inco, nbprop

	videagreg()

	// Wave de l'inconu de la mort qui tue
	Wave/Z dataindex,deltam, deltai, molecule, transindex, simu_proba, simu_mass, deltappm, test1, labelpro, incox,incoxdm,incoy, testlocal
	//Waves Text
	Wave/T list_prop, elements,List_formules, settings
	String lawave


	// SI il existe on détruit le panel3
	If(Wintype("panel3")!=0)
		Killwindow panel3
	endIf
	
	// La wave Progressbarre va contenir toutes les infos
	Make/O/N=(n,118)/D progressbarre
	progressbarre=1e12
	//execute "showprogbar()"
		for(k=0;k<n;k+=1)
			inco=roim[k]
			incox=roim[k]
			incoxdm=roidm[k]
			incoy=roii[k]
			analysemass("arg")
			Make/O/D/N=(nbprop) critere
			critere=1
				for(j=0;j<nbprop;j+=1)
					molecule=0
					tuelesspec()
					lawave="molprop_"+num2str(j)
					duplicate/O $lawave Molecule
					genestringmol()
					duplicate/O test1, testlocal
					testlocal= test1 * simu_proba[x] * abs( str2num(list_prop[j][1])- 1e6*deltam[x]/simu_mass[x])
					duplicate/O simu_proba transindex
					makeindex/R simu_proba transindex
						if(abs(labelpro[transindex[0]]-str2num(list_prop[j][1]))<0.1 && abs(str2num(list_prop[j][1]))<seuil && abs(labelpro[transindex[0]])<seuil)//0.1 arbitraire
							critere[j]=mean(testlocal)
						else
							critere[j]=1
						endif
				endfor
			// recup Proximite
			progressbarre[k][0]= wavemin(critere)
			// recup Masse inconnue
			progressbarre[k][1]= inco
			// recup de l'intensite de la masse inconnue
			progressbarre[k][2]= incoy[0]
			if(wavemin(critere)<1)
				Findvalue/V=(wavemin(critere)) critere
				molecule=0
				tuelesspec()
				lawave="molprop_"+num2str(v_value)
				duplicate/O $lawave Molecule
				// genere nom molecule
				genestringmol()
				//ajoutte affichage de masse trouvée et ajoute dans liste pour liste bag
				addthissimu()
				duplicate/O simu_proba transindex
				makeindex/R simu_proba transindex
				// recup le delta PPM
				progressbarre[k][3]= labelpro[transindex[0]]
				// recup des indices pour chaques éléments
				progressbarre[k][4,118]= molecule[y-4]
			else
				molecule=0
				tuelesspec()
				genestringmol()
				addthissimu()
				progressbarre[k][3]=0
				progressbarre[k][4,118]= molecule[y-4]
			endif
			doupdate
	endfor
	
	
	/////////////////////////////////////////////////////////
	// Sauvegarde du molBag correspondant à l'attribution
	// et des chargeBag, listBag, cibleBag associés.
	//
	// On place le résultat dans le molbag current
	// Possibilité d'imaginer un systeme sauvegarde
	// différent.
	String nom = "Last_Attribution"
	sauvegardeMolBagApresAttrib(nom)
	sauvegardeChargeBagApresAttrib(nom)
	sauvegardeListBagApresAttrib(nom)
	sauvegardeCibleBagApresAttrib(nom)
	videagreg()
	addthemolbag(nom)
	
	
	/////////////////////////////////////////////////////////
	// Fabrique la wave de Rapport
	//
	// Fabrique la wave somme sur les collones pour detecter les éléments
	// présent ou non et fait le tri
	Make/O/N=118 somme
	Matrixop/O somme=sumcols(progressbarre)
	for(k=0, j=0;k<118;k+=1)
		if(somme[0][j]==0)
			deletepoints/M=1 j,1, progressbarre, somme, titres
		else
			j+=1
		endif
	endfor
	// Fabrique les titres
	Make/O/T/N=(1,118) titres
	Titres[0][0]="Proximite"
	Titres[0][1]="MassInconnue"
	Titres[0][2]="Intensitepic"
	Titres[0][3]="Delta"
	Titres[0][4,117]=elements[y-4]
	// Fait le Rapport
	Make/O/T/N=(n+1,j+1) Rapport
	Rapport[0][0]="Formule"
	Rapport[0][1,j+1]=titres[0][y-1]
	Rapport[1,n][1,j+1]=num2str(progressbarre[x-1][y-1])
	Rapport[1,n][0]=List_formules[numpnts(List_formules)-n+x-1]
	string masseprecise
	for(k=1;k<n+1;k+=1)
		sprintf masseprecise "%.6f", progressbarre[k-1][1]
		Rapport[k][2]=masseprecise
	endfor
	
	/////////////////////////////////////////////////////////
	// Demarage compare Attrib
	//
	// Kill la windows de barre de progression
	If(Wintype("showprogbar")!=0)
		KillWindow showprogbar
	endIf
	
	//  KILL la window de comparaison d'attribution,
	// Ceci entraine un ménage dans les waves.
	If(Wintype("AttributionComparaison")!=0)
		KillWindow AttributionComparaison
	endIf
	
	// Lance les fonctions d'initialisation de compareAttrib
	prepareCompareAttribution(nom)
	removeFakes()
	// Démarre l'interface GUI de compare Attrib
	interfaceCompareAttribution()

	killwaves somme, titres
end //fullauto(nprochains,seuil)



Function GoForAttribROI(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar ppmthreshold
	switch( ba.eventCode )
		case 2: // mouse up
			fullautoroi(roi0,roi1,roi2,ppmthreshold)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


////////////////////// ORLEANS DATA LOADER   /////////////////////////////////



//serveur SFTP : sftp-lpc2e.cnrs-orleans.fr
//login : orbitrap
//pass : " "




function LoadProto()
	getchemin()
		
	GBLoadWave/Q/O/B/V/N=transient/T={2,2}/S=4/W=1/P=cheminfichier "signal"
	GBLoadWave/Q/O/B/V/N=size/T={4,2}/S=4/W=1 /P=cheminfichier "signal"
	wave size0
	wave transient0
	variable/G facteurreduc=1000
	variable/G facteurdigitalisation= size0(numpnts(size0)-1)/1e6
	DeletePoints numpnts(transient0)-2,2, transient0
	SetScale/P x 0,facteurdigitalisation,"s", transient0
	reduit(facteurreduc)
	creevisu()
end

function reduit(facteurreduc)
variable facteurreduc
variable/G facteurdigitalisation
wave transient0
	make/O/n=(numpnts(transient0)/facteurreduc) disptransient
	make/O/n=2 zonechoix
	disptransient=transient0[x*facteurreduc]
	zonechoix = transient0[x*facteurreduc]
	SetScale/P x 0,facteurdigitalisation*facteurreduc,"s", disptransient, zonechoix
end

function getchemin()
	pathinfo grandchemin
	NewPath /O/Q/Z cheminfichier 
endmacro



Function updatetransient(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	variable/G debuttransient
	variable/G fintransient
	variable/g hanningdo
	wave disptransient
	duplicate/O/R=((debuttransient/1000),(fintransient/1000)), disptransient, zonechoix,windowres,zonewindow
	if (hanningdo)
		mywindow(zonechoix)
	endif
	return 0
End

function mywindow(wavechoisie)
	wave wavechoisie	
	string/G nomfonction
	strswitch(nomfonction)	// numeric switch
		case "Hanning":
			WindowFunction Hanning,wavechoisie
			break
		case "Blackman":
			WindowFunction Blackman,wavechoisie
			break
		case "Cos1":
			WindowFunction Cos1,wavechoisie
			break
		case "cos2":
			WindowFunction cos2,wavechoisie
			break
		case "cos3":
			WindowFunction cos3,wavechoisie
			break
		case "cos4":
			WindowFunction cos4,wavechoisie
			break
		case "Hamming":
			WindowFunction Hamming,wavechoisie
			break
		case "Bartlett":
			WindowFunction Bartlett,wavechoisie
			break
		case "KaiserBessel20":
			WindowFunction KaiserBessel20,wavechoisie
			break
		case "KaiserBessel25":
			WindowFunction KaiserBessel25,wavechoisie
			break
		case "KaiserBessel30":
			WindowFunction KaiserBessel30,wavechoisie
			break
		case "Parzen":
			WindowFunction Parzen,wavechoisie
			break
		case "Poisson2":
			WindowFunction Poisson2,wavechoisie
			break
		case "Poisson3":
			WindowFunction Poisson3,wavechoisie
			break
		case "Poisson4":
			WindowFunction Poisson4,wavechoisie
			break
		case "Riemann":
			WindowFunction Riemann,wavechoisie
			break
		endswitch
end 

Function dohanning(name,value)
	String name
	Variable value
	variable/G facteurreduc
	variable/G debuttransient
	variable/G fintransient
	wave disptransient		
	strswitch (name)
		case "checkhanningon":
			mywindow(zonechoix)
			variable/G hanningdo=1
			CheckBox checkhanningof,value= 0
			CheckBox checkhanningon,value= 1
			showwindow()
			break
		case "checkhanningof":
			variable/G hanningdo=0
			duplicate/O/R=((debuttransient/1000),(fintransient/1000)), disptransient, zonechoix
			CheckBox checkhanningon,value= 0
			CheckBox checkhanningof,value= 1
			dowindow/K showwindo
			break
	endswitch
End

function showwindow() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/N=showwindo /W=(203,350,708,620) zonewindow
	ModifyGraph rgb(zonewindow)=(0,0,52224)
	wave disptransient
	GetAxis /W=visualyxan /Q bottom
	setaxis bottom V_min, V_max
	ControlBar 25
	selectwindow ("",1,"hanning")
	PopupMenu popup0,pos={20,0},size={247,21},title ="fenetre a utiliser"
	PopupMenu popup0,mode=1, proc=selectwindow
	PopupMenu popup0,value= #"\"Hanning;Blackman;Cos1;cos2;cos3;cos4;Hamming;Bartlett;KaiserBessel20;KaiserBessel25;KaiserBessel30;Parzen;Poisson2;Poisson3;Poisson4 ;Riemann \""
EndMacro

Function selectwindow (ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	string/G nomfonction= popStr
	wave zonewindow
	variable/G debuttransient
	variable/G fintransient
	duplicate/O/R=((debuttransient/1000),(fintransient/1000)), disptransient, zonechoix
	zonewindow=1
	execute "WindowFunction "+popStr +", zonewindow"
	execute "WindowFunction "+popStr +", zonechoix"	
End

function creevisu()
	PauseUpdate; Silent 1		
	dowindow/K visualyxan
	dowindow/K showwindo
	Display/N=visualyxan/W=(203,41,708,334) disptransient
	appendtograph zonechoix
	ModifyGraph axisEnab(left)={0,0.95}
	ModifyGraph rgb(zonechoix)=(0,0,52224)
	ControlBar 20
	Variable/G debuttransient, fintransient, hanningdo, zerofill
	SetVariable setdebut pos={0,2},size={110,20},title="debut (ms)"
	SetVariable setdebut, value=debuttransient,limits={0,2000,50} , proc=updatetransient
	SetVariable setdefin, pos={115,2},size={85,20},title="fin (ms)"
	SetVariable setdefin, value=fintransient,limits={0,2000,50}, proc=updatetransient
	SetVariable setzerofill, pos={200,2},size={70,20},title="Zer F"
	SetVariable setzerofill, value=zerofill,limits={0,5,1}
	CheckBox checkhanningon,pos={280,2},size={80,20},title="windowing ON",value= hanningdo,mode=0,proc=dohanning, side=1
	CheckBox checkhanningof,pos={345,2},size={50,20},title="OFF",value= (1-hanningdo),mode=0,proc=dohanning, side=1
	button doft, pos={400,0},size={75,20},title="DO FFT",proc=doFFT
	button nextfile, pos={560,0},size={75,20},title="Next",proc=loadnextalyxan
	button prevfile, pos={480,0},size={75,20},title="Prev",proc=loadprevalyxan
	pathinfo cheminfichier
	string info = "\\Z08"+S_path
	TextBox/C/N=text0/F=0/A=RC/X=-2/Y=51.5 info
	if (hanningdo)
		dohanning("checkhanningon",0)
	else
		dohanning("checkhanningof",1)
	endif
EndMacro

function loadnextalyxan(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2:
		dowindow/K 	visualyxan
		getcheminnext()
		
		GBLoadWave/Q/O/B/V/N=transient/T={2,2}/S=4/W=1/P=cheminfichier "signal"
		GBLoadWave/Q/O/B/V/N=size/T={4,2}/S=4/W=1 /P=cheminfichier "signal"
		wave size0
		wave transient0
		variable/G facteurreduc=1000
		variable/G facteurdigitalisation= size0(numpnts(size0)-1)/1e6
		DeletePoints numpnts(transient0)-2,2, transient0
		SetScale/P x 0,facteurdigitalisation,"s", transient0
		reduit(facteurreduc)
		creevisu()
		if (ba.eventMod ==8) 	// l'utilisateur a enfoncé la touche control	
			dofft("bobo")	
		endif
		break
	endswitch
end

function loadprevalyxan(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2:	
		dowindow/K 	visualyxan
		getcheminprev()
		
		GBLoadWave/Q/O/B/V/N=transient/T={2,2}/S=4/W=1/P=cheminfichier "signal"
		GBLoadWave/Q/O/B/V/N=size/T={4,2}/S=4/W=1 /P=cheminfichier "signal"
		wave size0
		wave transient0
		variable/G facteurreduc=1000
		variable/G facteurdigitalisation= size0(numpnts(size0)-1)/1e6
		DeletePoints numpnts(transient0)-2,2, transient0
		SetScale/P x 0,facteurdigitalisation,"s", transient0
		reduit(facteurreduc)
		creevisu()
		if (ba.eventMod ==8) 	// l'utilisateur a enfoncé la touche control	
			dofft("bobo")	
		endif
		break
	endswitch
end

function getcheminnext()
	pathinfo cheminfichier
	newpath /O/Q/Z upperfolder,  ParseFilePath(1, S_path, ":", 1, 0) 
	string listdetoutlesdossier =indexedDir (upperfolder, -1,0)
	pathinfo cheminfichier
	string actuel =ParseFilePath(0, S_path, ":", 1, 0) 
	variable indice =  WhichListItem(actuel,listdetoutlesdossier)  		
	pathinfo upperfolder
	string nextone=StringFromList((indice+1), listdetoutlesdossier)
	if (strlen(nextone)!=0)
		nextone= s_path+nextone
		NewPath /O/Q/Z cheminfichier nextone
	endif
end
function getcheminprev()
	pathinfo cheminfichier
	newpath /O/Q/Z upperfolder,  ParseFilePath(1, S_path, ":", 1, 0) 
	string listdetoutlesdossier =indexedDir (upperfolder, -1,0)
	pathinfo cheminfichier
	string actuel =ParseFilePath(0, S_path, ":", 1, 0) 
	variable indice =  WhichListItem(actuel,listdetoutlesdossier)  		
	pathinfo upperfolder
	string nextone = s_path+StringFromList((indice-1), listdetoutlesdossier)
	if (indice>=1)
		NewPath /O/Q/Z cheminfichier nextone
	endif
end

function doclean() : ButtonControl
	variable/G facteurreduc
	variable/G debuttransient
	variable/G fintransient
	variable/G facteurdigitalisation
	Variable/G masseexacte
	variable/G frequenceexacte
	variable/G hanningdo
	wave ftmasse
	wave ftintensite
	findlevel/P/Q ftmasse 500
	deletepoints 0,(v_levelX), ftmasse,ftintensite
	variable minimum=wavemin(ftmasse)
	variable maximum,pointnimum,pointmaximum
	do
		findlevel/P/Q ftmasse minimum
		pointnimum=v_levelX
		maximum = 0.9987*floor(minimum+1)
		findlevel/P/Q ftmasse maximum
		pointmaximum=v_levelX
		if (pointmaximum<pointnimum)
			deletepoints pointmaximum,(pointnimum-pointmaximum), ftmasse,ftintensite
			minimum=1.0015*ceil (maximum)
		else
			break
		endif
	while (1)
end

function dofft(ctrlName) : ButtonControl
	String ctrlName
	variable/G facteurreduc
	variable/G debuttransient
	variable/G fintransient
	variable/G facteurdigitalisation
	Variable/G masseexacte
	variable/G frequenceexacte
	variable/G hanningdo
	variable/G zerofill
	wave transient0
	duplicate/S/O/R=((debuttransient/1000),(fintransient/1000)) transient0 wavesourse,windowres
	if (mod(numpnts(wavesourse), 2)!=0)
		DeletePoints numpnts(wavesourse)-1,2, wavesourse
	endif
	if (hanningdo)
		mywindow(wavesourse)
	endif
	insertpoints 0,zerofill*(numpnts(wavesourse)), wavesourse
	SetScale/P x 0,facteurdigitalisation,"s",wavesourse
	
	FFT/OUT=3/DEST=FTMS wavesourse
	duplicate/O FTMS FTmasse FTintensite
	
	FTmasse=masseexacte*(frequenceexacte/x)^2
	showfftresult()
end

function showfftresult() : Graph
	PauseUpdate; Silent 1		// building window...
	dowindow/F FTmas
	if (v_flag!=1)
		Display /N=FTmas/W=(714,351,1224,620) FTintensite vs FTmasse
		ModifyGraph log(left)=1
		SetAxis bottom 0,300
		ModifyGraph axisEnab(left)={0,0.9}
//		button docle, pos={400,0},size={100,20},title="clean chem",proc=dodoclean
		button doscal, pos={200,0},size={100,20},title="scale 0-500",proc=domyscal
		button doaddelem, pos={20,0},size={165,20},title="eval element",proc=doaddelementonFFT
		button goattrb, pos={400,0},size={100,20},title="go attributor",proc=goattrib
		
	endif
	dowindow/F FTresult
	if (v_flag!=1)
		Display/N=FTresult /W=(714,41,1224,332) FTMS
		ModifyGraph log(left)=1
		Cursor A, FTMS, wavemax(FTMS)
		ControlBar 20
		Variable/G masseexacte
		variable/G frequenceexacte
		SetVariable masex pos={0,2},size={200,20},title="masse exacte"
		SetVariable masex, value=masseexacte , proc=updatemasscale
		SetVariable frequenceex, pos={200,2},size={200,20},title="fréquence exacte"
		SetVariable frequenceex, value=frequenceexacte, proc=updatemasscale
		button goattrib, pos={420,0},size={100,20},title="calibre sur élément",proc=usesimu
		
	endif
End

function domyscal(ctrlName) : ButtonControl
	String ctrlName
	dowindow/F FTmas
	SetAxis bottom 0,500
end

function doaddelementonFFT(ctrlName) : ButtonControl
	String ctrlName
	wave FTmasse, FTintensite, transient0,size0,disptransient,zonechoix,wavesourse,FTMS, simu_mass, simu_proba
	SetScale/P x 0,1,"", FTintensite,FTmasse	
	dowindow/F FTmas
	wavestats/Q/R =[numpnts(FTintensite)*0.8,numpnts(FTintensite)] FTintensite
	//wavestats/Q  FTintensite
	variable noiselevel=V_avg
	genestringmol()
	//print noiselevel
	wavestats/Q simu_proba
	findlevel/Q FTmasse, simu_mass[v_maxloc]
	variable maxloc=FTintensite[V_levelx]
	variable intref=simu_proba[v_maxloc]
	simu_proba*=(maxloc-noiselevel)/intref
	simu_proba+=noiselevel
	if(strlen( WaveList("simu_proba", ";","WIN:FTmas"))==0)
		appendtograph simu_proba vs simu_mass
	endif
	ModifyGraph mode(simu_proba)=3,rgb(simu_proba)=(0,0,52224)
	SetAxis bottom (wavemin(simu_mass)-1), (wavemax(simu_mass)+1)
	SetAxis left noiselevel, (wavemax(simu_proba)*2)

end

function usesimu(ctrlName) : ButtonControl
	String ctrlName
	wave FTmasse, FTintensite, transient0,size0,disptransient,zonechoix,wavesourse,FTMS, simu_mass, simu_proba
	variable masseref
	wavestats/Q simu_proba
	masseref= simu_mass(v_maxloc)
	variable freqref
	freqref=xcsr(A)
	Variable/G masseexacte, frequenceexacte
	masseexacte=masseref
	frequenceexacte=freqref
	FTmasse=masseexacte*(frequenceexacte/x)^2
	dowindow/F FTmas
	SetAxis bottom (wavemin(simu_mass)-1), (wavemax(simu_mass)+1)
	if(strlen( WaveList("simu_proba", ";","WIN:FTmas"))==0)
		appendtograph simu_proba vs simu_mass
		ModifyGraph mode(simu_proba)=3,rgb(simu_proba)=(0,0,52224)
	endif
end

function goattrib(ctrlName) : ButtonControl
	String ctrlName
	wave FTmasse, FTintensite, transient0,size0,disptransient,zonechoix,wavesourse,FTMS,zonewindow
	variable/G zerofill 
	SetScale/P x 0,1,"", FTintensite,FTmasse
	smooth 5,FTintensite 
	Resample/DOWN=(zerofill+1) FTintensite,FTmasse
	SetScale/P x 0,1,"", FTintensite,FTmasse
	doclean()
	sort ftintensite,ftmasse, ftintensite
	wavestats/Q ftintensite
	findlevel/Q ftintensite, (V_avg*0.8)
	deletepoints 0, V_levelX, ftintensite, ftmasse
	sort  ftmasse, ftintensite, ftmasse
	duplicate/O ftintensite CosmOrbi_1stripped
	duplicate/O ftmasse CosmOrbi_0stripped
	listelesdonnees()
	dowindow/F advancedmanager
	dowindow/B visualyxan
	dowindow/B FTresult
	dowindow/B FTmas
	dowindow/K showwindo
	//deletenoiseFAT(FTintensite,FTmasse)
	//deletenoiseFINE(FTintensite,FTmasse)
	//dowindow/K visualyxan
	//dowindow/K FTresult
	//dowindow/K FTmas
	//dowindow/K showwindo
	//KillWaves  transient0,size0,disptransient,zonechoix,wavesourse,FTMS,FTintensite, FTmasse,zonewindow 
end

function dodoclean(ctrlName) : ButtonControl
	String ctrlName
 	doclean()
end

Function updatemasscale(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	Variable/G masseexacte, frequenceexacte
	wave FTmasse
	FTmasse=masseexacte*(frequenceexacte/x)^2
	return 0
End



Window laura() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(714,350.75,1224,620) FTintensite vs FTmasse
	ModifyGraph log(left)=1
	SetAxis bottom 0,500
EndMacro

============
//modif janvier 2015 par Roland pour downloader les donnees depuis le FTP de l'IPAg //Revue le 10 juillet 2015 
function loadIPAGdata() 
string url="ftp://janus.obs.ujf-grenoble.fr/textified"
string nom="04029b"
string urlDAT
String fileName="info.txt"
//browseURLURL
Prompt nom,"Enterdataname:"
doprompt "indicatethefiletodownload",nom
if(V_Flag)
     return 0//usercanceled
endif
url+="/cahier"+nom[1]+"/"+nom
urlDAT=url+"/info.txt"
//printurldat
FTPDownload/P=igor/O=1/T=1/V=2/U="orbitrap"/W="titan"urlDAT,"info.txt"
if (v_flag)
     return 0
endif
print S_Filename
loadwave/P=igor/Q/O/N=infodataIPAG/J/D/K=2"info.txt"
getdataIPAG(url,nom)
deletefile/Z/P=igor"info.txt"
deletefile/Z/P=igor"data.txt"
listelesdonnees()
end


//modif janvier 2015 par Roland pour downloader les donnees depuis le FTP de l'IPAg //Revue le 10 juillet 2015  
function getdataIPAG(url,nom) 
string url,nom 
wave/T wav=infodataipag0 
variable nombrefiltre 
variable i 
string lis_filtre=""
variable datanumber
nombrefiltre=0.5*(numpnts(wav)-1)
for(i=0;i<nombrefiltre;i+=1)
     lis_filtre+=wav[1+2*i]
     lis_filtre+=";"
endfor
//printlis_filtre
Prompt datanumber,"filtre",popup,lis_filtre
DoPrompt "choosethefiltertodownload",datanumber
if(V_Flag)
return 0//usercanceled
endif
string urldat=url+"/dataFatNo"+num2str(datanumber)+".txt"
string fichiery,fichierx
FTPDownload/P=igor/O=1/T=1/V=2/U="orbitrap"/W="titan"urlDAT,"data.txt"
LoadWave/P=igor/O/A/Q/G/D"data.txt"
fichiery=nom+num2str(datanumber)+"_1stripped"
fichierx=nom+num2str(datanumber)+"_0stripped"
Rename wave0,$fichierx
Rename wave1,$fichiery
end


Menu "Load Waves"
                "Load IPAG file", loadIPAGdata()
                "Load Cosmorbitrap file", LoadProto()
                submenu "COSIMA"
                "Load COSIMA single file", loadcosima()
                "Load all COSIMA TABS from a folder" , massloadCosimaTABS()
                end
End

//========MSMS ROLAND STARTS=============
Function gotomsms(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			modifysetupformsm()
			break
	endswitch

	return 0
End

function modifysetupformsm()
	string existe
	existe=winlist("panel10",",", "")
	if (cmpstr(existe,"")==0)
		paneldix()
	else
		dowindow/F panel10
	endif
	existe=winlist("loadedData",",", "")
	if (strlen(existe)==0)
		nb1()
	else
		dowindow/F loadedData
	endif
	existe=winlist("evaluateur",",", "")
	if (strlen(existe)==0)
		execute "evaluateur()"
	else
		dowindow/F evaluateur
	endif
	dowindow/F mendeleiev
	Button msmser,pos={833,157},size={102,25},proc=backfrommsms,title="backtoattributor"
	return 0
end

function nb1() : Panel
	PauseUpdate; Silent 1		// building window...
	String nb = "LoadedData"
	NewNotebook/N=$nb/F=1/V=1/W=(830,37,1237,720)
	
	Notebook $nb defaultTab=20, statusWidth=238, backRGB=(57344,65280,48896), pageMargins={71,71,71,71}
	Notebook $nb showRuler=0, rulerUnits=2, updating={1, 60}
	Notebook $nb newRuler=Normal, justification=0, margins={0,0,468}, spacing={0,0,0}, tabs={}, rulerDefaults={"Arial",12,1,(0,0,0)}
	Notebook $nb newRuler=title, justification=0, margins={0,0,468}, spacing={0,0,0}, tabs={}, rulerDefaults={"Arial",10,3,(0,0,0)}
	Notebook $nb newRuler=insert, justification=0, margins={15,15,468}, spacing={0,0,0}, tabs={}, rulerDefaults={"Arial",10,0,(0,26112,0)}
	Notebook $nb ruler=Normal, text="MS/MS spectra\r\r"
	return 0
End

Window evaluateur() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(0,37,827,579)
	make/O/N=2 lamassemsms lintensite, massesimu, fakeint
	append lintensite vs lamassemsms
	append/L=simu fakeint vs massesimu
	ModifyGraph noLabel(simu)=2,axisEnab(left)={0.05,1},axisEnab(simu)={0,0.05}
	ModifyGraph freePos(simu)=0
	ModifyGraph mode(fakeint)=1,rgb(fakeint)=(0,0,0)
	SetAxis simu 0,1
	lintensite=1
	fakeint=1
	controlbar 28
	SetVariable setcharge,pos={1,2},size={60,26},proc=Setatominmolec,title="+/-",limits={-1,1,2 }, fsize=12
	SetVariable setcarbone,pos={70,2},size={60,26},proc=Setatominmolec,title="C",limits={0,50,1 },fsize=12
	SetVariable sethydrogen,pos={140,2},size={60,26},proc=Setatominmolec,title="H",limits={0,80,1 },fsize=12
	SetVariable setnitrogen,pos={210,2},size={60,26},proc=Setatominmolec,title="N",limits={0,45,1 },fsize=12
	SetVariable setoxygen,pos={280,2},size={60,26},proc=Setatominmolec,title="O",limits={0,45,1 },fsize=12
	Button addDer,pos={350,2},size={60,26},proc=addDeuterium,title="With S?"
	Button addNaer,pos={420,2},size={60,26},proc=addsodium,title="With Na?"
	valdisplay deltaenppm,pos={490,2},size={200,26},title="delta ppm", value= K14, fsize=12
	SetVariable setparent,pos={700,2},size={200,26},title="MassParent",fsize=12, disable=1
	SetVariable setcharge,value= K9
	SetVariable setcarbone,value= K10
	SetVariable sethydrogen,value= K11
	SetVariable setnitrogen,value= K12
	SetVariable setoxygen,value= K13
	SetVariable setparent,value= K15
	variable/G withD=0
EndMacro

Function addDeuterium(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			KillControl addDer
			SetVariable setdeuterium,pos={350,2},size={60,26},proc=Setatominmolec,title="S",limits={-1,45,1 },fsize=12
			SetVariable setdeuterium,value= K19
			variable/G withD=1
			break
	endswitch
	return 0
End
 
Function backfromaddDeuterium() : ButtonControl
	KillControl setdeuterium
	Button addDer,pos={350,2},size={60,26},proc=addDeuterium,title="With D?"
	variable/G withD=0
	return 0
End

Function addsodium(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			KillControl addNaer
			SetVariable setSodium,pos={420,2},size={60,26},proc=Setatominmolec,title="Na",limits={-1,1,2 },fsize=12
			SetVariable setSodium,value= K18
			variable/G withNa=1
			break
	endswitch
	return 0
End
 
Function backfromaddsodium() : ButtonControl
	KillControl setSodium
	Button addNaer,pos={420,2},size={60,26},proc=addsodium,title="With Na?"
	variable/G withNa=0
	return 0
End


function Setatominmolec(ctrlName,varNum,varStr,varName) :  SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
	Make/O/D/N=0 massesimu, fakeint, nombC, nombH, nombN, nombO,nombD,nombNa
	wave/Z mendev_masses
	variable setcharge=K9
	variable Nc=K10, indiceC=0
	variable NH=K11, indiceH=0
	variable Nn=K12, indiceN=0
	variable No=K13, indiceO=0
	variable ND=K19, indiceD=0
	variable NNa=K18, indiceNa=1
	variable DOUpa 
	nvar withD
	if (withD==0)
		 ND=0
	endif
	if (ND==-1)
		 backfromaddDeuterium()
		 K19=0
		 return 0
	endif
	nvar withNa
	if (withNa==0)
		 NNa=0
		 indiceNa=0
	endif
	if (NNa==-1)
		 backfromaddsodium()
		 K18=0	 
		 return 0
	endif
	do
		do
			do
				do
					do
						insertpoints 0,1, massesimu, fakeint, nombC, nombH, nombN, nombO,nombD,nombNa
						//massesimu[0]=-0.0005*setcharge+indiceC*12+indiceH*mendev_masses[0][0]+indiceD*mendev_masses[0][1]+indiceN*mendev_masses[6][0]+indiceO*mendev_masses[7][0]+indiceNa*mendev_masses[10][0]
						massesimu[0]=-0.0005*setcharge+indiceC*12+indiceH*mendev_masses[0][0]+indiceD*mendev_masses[15][0]+indiceN*mendev_masses[6][0]+indiceO*mendev_masses[7][0]+indiceNa*mendev_masses[10][0]	// janvier 2016 pour utiliser la masse du soufre au lieu de celle du deuterium  
						 nombC[0]=indiceC
						 nombH[0]=indiceH
						 nombN[0]=indiceN
						 nombO[0]=indiceO
						 nombD[0]=indiceD
						 nombNa[0]=indiceNa
						 indiceD+=1
					while (indiceD<=Nd)
					indiced=0
				indiceH+=1
				while (indiceH<=NH)
				indiceH=0
				indiceC+=1
			while (indiceC<=Nc)
			indiceC=0
			indiceN+=1
		while (indiceN<=Nn)
		indiceN=0
		indiceO+=1
	while (indiceO<=No)
	indiceO=0
	sort/R massesimu massesimu, fakeint, nombC, nombH,nombD, nombN, nombO,nombNa
	fakeint=1	
	if (k15==0)
		wave/Z w_coef
		K14=1e6*(w_coef[2]-massesimu[0])/massesimu[0]
	else
		K14=1e6*(k15-massesimu[0])/massesimu[0]
	endif
	if (varnum!=111)
		execute "zoom()"
	endif
end


Function backfrommsms(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			modifysetupfbackfrommsm()
			break
	endswitch

	return 0
End

function modifysetupfbackfrommsm()
	wave/Z lamassemsms
	wave/Z lintensite
	string fichiery
	string fichierx
	dowindow/B evaluateur
	dowindow/B panel10
	dowindow/B Loadeddata
	Button msmser,pos={833,157},size={102,25},proc=gotoMSMS,title="doMSMS"
	if (numpnts(lintensite)>2)
		fichiery="msmsSpectrum_1stripped"
		fichierx="msmsSpectrum_0stripped"
		duplicate/O lintensite $fichiery
		duplicate/O lamassemsms $fichierx
		listelesdonnees()
	endif
	return 0
end
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function paneldix() : Panel			//contiendra les données chargée
	PauseUpdate; Silent 1		// building window...
	NewPanel /N=panel10/W=(954,797,1103,979)/N=panel10
	ModifyPanel cbRGB=(0,26112,0)
	Button chargeTout,pos={10,5},size={80,20},proc=loadmsms,title="load data"
	Button erasebutton,pos={95,5},size={45,20},proc=killallMSMSdata,title="erase"
	Button MontreTout,pos={10,30},size={80,20},proc=montredata,title="show data"
	Button logit,pos={95,30},size={45,20},proc=gotolog,title="log"
	Button usefull,pos={10,55},size={80,20},proc=cleanup,title="clean HCNO"
	Button noparentbutton,pos={95,55},size={45,20},proc=setparent,title="Set P"
	Button zom,pos={10,80},size={80,20},proc=zoomzoom,title="zoom it"
	Button calespec,pos={95,80},size={45,20},proc=docalib,title="recale"
	Button evalue,pos={10,105},size={130,20},proc=doevaluetaion,title="eval MSMS"
	Button gauchebutton,pos={10,130},size={15,20},proc=zoomgauche,title="<"
	Button checkbutton,pos={30,130},size={40,20},proc=checkmsms,title="check"
	Button droitebutton,pos={75,130},size={15,20},proc=zoomdroit,title=">"
	Button killbutton,pos={95,130},size={45,20},proc=killresult,title="kill"
	Button savemsmsbutton,pos={10,155},size={130,20},proc=savemsms,title="save MSMS"
End

Function loadmsms(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			execute "faitledossier()"
			break
	endswitch
	return 0
End

macro faitledossier()
silent 1
	variable indice
	string cheminfersfichier
	GetFileFolderInfo/Q
	NewPath/O/Q leGRANDdossier ParseFilePath(1, S_Path, ":", 1,1)
	do
		cheminfersfichier=IndexedDir(leGRANDdossier, indice, 0)
		chargelesdatatext(cheminfersfichier)
		indice+=1
		cheminfersfichier=IndexedDir(leGRANDdossier, indice, 0)
	while (strlen(cheminfersfichier)!=0)
	dowindow/F panel10
endmacro

macro chargelesdatatext(cheminversfichier)
	string cheminversfichier
	silent 1
	make/O/N=2/T typems=""
	make/O/N=2/T debutfiltre=""
	variable nombrefiltre 
	variable indicefiltre = 0
	string nomfichierdata, nomint, nommass
	if (strlen(cheminversfichier)==0)
		GetFileFolderInfo/Q
		NewPath/O/Q ledossier ParseFilePath(1, S_Path, ":", 1, 0)
	else
		pathinfo leGRANDdossier
		S_Path+=cheminversfichier
		NewPath/O/Q ledossier s_path
	endif
	debutfiltre(0)=cheminversfichier
	debutfiltre(1)=S_Path
	LoadWave/Q/N=Info/P=ledossier/J/O/K=2 "info.txt"
	//nombrefiltre = str2num (info0(0))
	nombrefiltre = (numpnts (info0)-1)/2
	Notebook  LoadedData ruler=title, text="\r#"+cheminversfichier+"\r"
	do
		nomfichierdata="data"+num2str (indicefiltre+1)+".txt"
		nomint="int"+cheminversfichier+num2str(indicefiltre+1)
		nommass="masse"+cheminversfichier+num2str(indicefiltre+1)
		LoadWave/Q/J/D/N=dat/P=ledossier nomfichierdata
		duplicate/O dat1,$nomint
		duplicate/O dat0,$nommass; 
		//Notebook  LoadedData ruler=insert, text=num2str(indicefiltre+1)+"_"+info0(indicefiltre*2+2)+"\r"
		Notebook  LoadedData ruler=insert, text=num2str(indicefiltre+1)+"_"+info0(indicefiltre*2+1)+"\r"
		indicefiltre+=1
	while (indicefiltre<nombrefiltre)
	killwaves dat0, dat1
endmacro

Function killallMSMSdata(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			string laliste1=wavelist("int0*",";", "")
			//print laliste
			string laliste2=wavelist("masse0*",";", "")
			variable i
			string tracename, tracename2
			do
				traceName= StringFromList(i,laliste1)
				traceName2= StringFromList(i,laliste2)
				if( strlen(traceName) == 0 )
					break
				endif
				killwaves $tracename $tracename2
				i += 1
			while (1)	// exit is via break statement
			dowindow/K loadedData
			nb1()
			break
	endswitch
	return 0
End


Function montredata(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			execute "montremsms()"
			break
	endswitch
	return 0
End

macro montremsms()
	variable maximum, indice
	string nomtag
	variable/G testzoom=0
	dowindow/F evaluateur
	SetVariable setparent disable=1
	K15=0
	if (strlen(WaveList("resultI", "", "WIN:evaluateur" ))>0)
		remove resultI
		do
			nomtag="text"+num2str(indice)
			TextBox/K/N=$nomtag
			indice+=1
		while (indice<(numpnts(resultI)))
	endif
	if (strlen(WaveList("resultI", "", "WIN:evaluateur" ))>0)
		RemoveFromGraph resultI
		RemoveFromGraph lintensite#1
		RemoveFromGraph fakeint#1
	endif
	ModifyGraph axisEnab(bottom)={0,1}
	GetSelection notebook,LoadedData,2
	analysenomfichier (S_selection)
	dowindow/F evaluateur
	duplicate/O lintensite intensite_Hist logintensite
	logintensite=log(lintensite+1)
	Histogram/B=4 logintensite,intensite_Hist
	CurveFit/Q/NTHR=1 gauss  intensite_Hist
	variable bruitmoyen=10^(w_coef(2)), largeur= (10^(w_coef(2)+w_coef(3)/2)-bruitmoyen)
	DrawAction delete
	SetDrawEnv ycoord= left//, dash= 1;DelayUpdate
	DrawLine 0,bruitmoyen,1,bruitmoyen
	SetDrawEnv ycoord= left;DelayUpdate
	DrawLine 0,(bruitmoyen+2*largeur),1,(bruitmoyen+2*largeur)
	wavestats/Q lintensite
	maximum = V_max
	if ((bruitmoyen-largeur)>1)
		SetAxis left bruitmoyen-largeur,maximum
	else
		SetAxis/A left
	endif
	wavestats/Q lamassemsms
	SetAxis bottom V_min,V_max
	dowindow/F panel10
end

macro analysenomfichier(selectionotebook)
	string selectionotebook, nomfichier, nomintensite, nommasse,cheminsauvegarde
	variable indice,test
	string/G chiffre=""
	do
		chiffre+=selectionotebook[indice]
		indice+=1
		test=cmpstr(selectionotebook[indice+1],"_")
	while (test==0)
	notebook LoadedData findText={"#", 16 }
	notebook LoadedData selection={startOfParagraph, endOfChars }
	GetSelection notebook,LoadedData,2
	nomfichier=S_selection[1, 1000]
	notebook LoadedData findText={chiffre+"_", 0 }
	notebook LoadedData selection={startOfParagraph, endOfChars }
	pauseupdate
	nomintensite="int"+nomfichier+chiffre
	nommasse="masse"+nomfichier+chiffre
	duplicate/O $nomintensite lintensite
	duplicate/O $ nommasse lamassemsms
	lintensite+=1
	sort lintensite lintensite lamassemsms
	findlevel/Q  lintensite, 1.1
	deletepoints 0,v_levelX, lintensite, lamassemsms
	sort lamassemsms lamassemsms lintensite
	pathinfo leGRANDdossier
	cheminsauvegarde=s_path
	cheminsauvegarde+=nomfichier
	NewPath/O/Q lasauvegarde cheminsauvegarde
end

Function gotolog(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			dowindow/F evaluateur
			ModifyGraph log(left)=1
			dowindow/F panel10
			Button logit,proc=gotolin,title="lin"
			dowindow/F panel10
			break
	endswitch
	return 0
End

Function gotolin(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			dowindow/F evaluateur
			ModifyGraph log(left)=0
			dowindow/F panel10
			Button logit,proc=gotolog,title="log"
			dowindow/F panel10
			break
	endswitch
	return 0
End

Function cleanup(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			execute "usefulmass()"
			break
	endswitch
	return 0
End

macro usefulmass()
	variable indicemasse
	variable limitebasse limitehaute, pointbas, pointhaut,finspectre
	silent 1
	pauseupdate
	indicemasse=lamassemsms(0)
	finspectre=lamassemsms(numpnts(lamassemsms)-1)
	do
		//print indicemasse
		findlevel/Q lamassemsms, (indicemasse*1.0016)
		pointbas= v_levelX
		findlevel/Q lamassemsms, (1+indicemasse* 0.9996)
		pointhaut= v_levelX
		deletepoints pointbas, (pointhaut-pointbas), lintensite
		deletepoints pointbas, (pointhaut-pointbas), lamassemsms
		indicemasse+=1
	while (indicemasse<finspectre)
	dowindow/F panel10
endmacro

Function setparent(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z massesimu
			dowindow/F evaluateur
			SetVariable setparent disable=0
			K15=massesimu[0]
			break
	endswitch

	return 0
End

Function zoomzoom(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	nvar testzoom
	switch( ba.eventCode )
		case 2: // mouse up
			if (testzoom==0)
				execute "zoom()"
				else
				if (testzoom==1)
					execute "zoomparent()"
				else
					execute "zoomparentstrong()"
				endif
			endif
			break
	endswitch
	return 0
End

macro zoom ()
	GetSelection notebook,LoadedData,2
	analysestring (S_selection)
	dowindow/F evaluateur
	wavestats/Q lamassemsms
	if (k15==0)
		setaxis bottom V_min, (str2num(nomion)+5)
	else
		setaxis bottom V_min, (k15+5)
	endif
	dowindow/F panel10
	testzoom=1
end

macro zoomparent ()
	GetSelection notebook,LoadedData,2
	analysestring (S_selection)
	dowindow/F evaluateur
	if (k15==0)
		setaxis bottom (str2num(nomion)-2),(str2num(nomion)+2)
	else
		setaxis bottom (floor(k15)-2),(floor(k15)+2)
	endif
	dowindow/F panel10
	testzoom=2
end

macro zoomparentstrong ()
	variable pointhaut pointbas
	dowindow/F evaluateur
	if (k15==0)
		GetSelection notebook,LoadedData,2
		analysestring (S_selection)
		findlevel/Q lamassemsms, (str2num(nomion)-0.2)
		pointbas= v_levelX
		findlevel/Q lamassemsms, (str2num(nomion)+0.2)
		pointhaut= v_levelX
		CurveFit/Q/NTHR=0 gauss  lintensite[pointbas,pointhaut] /X=lamassemsms 
		setaxis bottom (W_coef[2]*0.9999),(W_coef[2]*1.0001)
	else
		setaxis bottom (k15*0.9999),(k15*1.0001)
	endif
	Setatominmolec("",111,"","")
	dowindow/F panel10
	testzoom=0
end

macro analysestring(filtre)
	string filtre, caractere,  
	variable indice, indicestring,positionparenthese, indiceMS, premier, parent, fille
	indice+=1
	silent 1
	do
		indicestring+=1
		caractere= filtre[indicestring]
	while (cmpstr (caractere, "[") !=0)
	positionparenthese = indicestring
	do
		indicestring-=1
		caractere= filtre[indicestring]
		if (indicestring==0)
			break
		endif
	while (cmpstr (caractere, "@") !=0)
	string/G nomion = filtre[(indicestring-6),(indicestring-1)]
	string/G energie = filtre[(indicestring+4),positionparenthese-2]
endmacro

Function docalib(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			 execute "zoomparentstrong ()"
			wave/Z lamassemsms
			wave/Z massesimu
			wave/Z w_coef
			variable massetheorique
			variable massemesuree
			massetheorique=massesimu[0]
			if (K15==0)
				massemesuree=W_coef[2]
			else
				massemesuree=k15
			endif
			lamassemsms+=(massetheorique-massemesuree)
			Setatominmolec("",111,"","")
			dowindow/F panel10
			break
	endswitch
	return 0
End

Function doevaluetaion(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			execute "evalMSMS ()"
			break
	endswitch
	return 0
End

macro evalMSMS ()
	silent 1
	dowindow/F evaluateur
	if (strlen(WaveList("resultI", ";", "WIN:evaluateur"))!=0)
		remove resultI
	endif
	variable v_fiterror v_fitoptions indice, massemin, indiceresult pointhaut pointbas
	string nomtag, letag
	CurveFit/Q/NTHR=1 gauss  intensite_Hist
	variable bruitmoyen=10^(w_coef(2)), largeur= (10^(w_coef(2)+w_coef(3)/2)-bruitmoyen)
	if ((withD==0)&&(withNa==0))
		make/O/N=(2,8) resultmsms
		massemin= lamassemsms(0)
		do
			findlevel/Q lamassemsms, (massesimu(indice)*0.99995)
			pointbas= v_levelX
			findlevel/Q lamassemsms, (massesimu(indice)*1.00005)
			pointhaut= v_levelX
			if ((pointhaut-pointbas)>5)
				v_fitoptions=4
				v_fiterror=0
				CurveFit/Q/NTHR=0 gauss  lintensite[pointbas,pointhaut] /X=lamassemsms 
			 	if (V_fiterror==0)
			 		if (W_coef[1]>0)
			 			if (abs(W_coef[1]-W_coef[0])>(bruitmoyen+2*largeur))
							if ((abs(W_coef[2]-massesimu(indice))/massesimu(indice))<5e-6)
								insertpoints 0,1, resultmsms
								resultmsms[0][0]=W_coef[1]-W_coef[0]
								resultmsms[0][1]=W_coef[2]
								resultmsms[0][2]=nombC[indice]
								resultmsms[0][3]=nombH[indice]
								resultmsms[0][4]=nombN[indice]
								resultmsms[0][5]=nombO[indice]
								resultmsms[0][6]=W_sigma[1]
								resultmsms[0][7]=W_sigma[2]
								indiceresult+=1
							endif
						endif
					endif
				endif
			endif
			indice+=1
		while (massesimu(indice)>massemin)
	else
		if (withNa==0)
			make/O/N=(2,9) resultmsms
			massemin= lamassemsms(0)
			do
				findlevel/Q lamassemsms, (massesimu(indice)*0.99995)
				pointbas= v_levelX
				findlevel/Q lamassemsms, (massesimu(indice)*1.00005)
				pointhaut= v_levelX
				if ((pointhaut-pointbas)>5)
					v_fitoptions=4
					v_fiterror=0
					CurveFit/Q/NTHR=0 gauss  lintensite[pointbas,pointhaut] /X=lamassemsms 
					 if (V_fiterror==0)
					 	if (W_coef[1]>0)
					 		if (abs(W_coef[1]-W_coef[0])>(bruitmoyen+2*largeur))
								if ((abs(W_coef[2]-massesimu(indice))/massesimu(indice))<5e-6)
									insertpoints 0,1, resultmsms
									resultmsms[0][0]=W_coef[1]-W_coef[0]
									resultmsms[0][1]=W_coef[2]
									resultmsms[0][2]=nombC[indice]
									resultmsms[0][3]=nombH[indice]
									resultmsms[0][4]=nombD[indice]
									resultmsms[0][5]=nombN[indice]
									resultmsms[0][6]=nombO[indice]
									resultmsms[0][7]=W_sigma[1]
									resultmsms[0][8]=W_sigma[2]
									indiceresult+=1
								endif
							endif
						endif
					endif
				endif
				indice+=1
			while (massesimu(indice)>massemin)
		endif
		if (withD==0)
			make/O/N=(2,9) resultmsms
			massemin= lamassemsms(0)
			do
				findlevel/Q lamassemsms, (massesimu(indice)*0.99995)
				pointbas= v_levelX
				findlevel/Q lamassemsms, (massesimu(indice)*1.00005)
				pointhaut= v_levelX
				if ((pointhaut-pointbas)>5)
					v_fitoptions=4
					v_fiterror=0
					CurveFit/Q/NTHR=0 gauss  lintensite[pointbas,pointhaut] /X=lamassemsms 
					 if (V_fiterror==0)
				 		if (W_coef[1]>0)
				 			if (abs(W_coef[1]-W_coef[0])>(bruitmoyen+2*largeur))
								if ((abs(W_coef[2]-massesimu(indice))/massesimu(indice))<5e-6)
									insertpoints 0,1, resultmsms
									resultmsms[0][0]=W_coef[1]-W_coef[0]
									resultmsms[0][1]=W_coef[2]
									resultmsms[0][2]=nombC[indice]
									resultmsms[0][3]=nombH[indice]
									resultmsms[0][4]=nombN[indice]
									resultmsms[0][5]=nombO[indice]
									resultmsms[0][6]=nombNa[indice]
									resultmsms[0][7]=W_sigma[1]
									resultmsms[0][8]=W_sigma[2]
									indiceresult+=1
								endif
							endif
						endif
					endif
				endif
				indice+=1
			while (massesimu(indice)>massemin)
		endif
		if ((withD!=0)&&(withNa!=0))

			make/O/N=(2,10) resultmsms
			massemin= lamassemsms(0)
			do
				findlevel/Q lamassemsms, (massesimu(indice)*0.99995)
				pointbas= v_levelX
				findlevel/Q lamassemsms, (massesimu(indice)*1.00005)
				pointhaut= v_levelX
				if ((pointhaut-pointbas)>5)
					v_fitoptions=4
					v_fiterror=0
					CurveFit/Q/NTHR=0 gauss  lintensite[pointbas,pointhaut] /X=lamassemsms 
					if (V_fiterror==0)
				 		if (W_coef[1]>0)
				 			if (abs(W_coef[1]-W_coef[0])>(bruitmoyen+2*largeur))
								if ((abs(W_coef[2]-massesimu(indice))/massesimu(indice))<5e-6)
									insertpoints 0,1, resultmsms
									resultmsms[0][0]=W_coef[1]-W_coef[0]
									resultmsms[0][1]=W_coef[2]
									resultmsms[0][2]=nombC[indice]
									resultmsms[0][3]=nombH[indice]
									resultmsms[0][4]=nombD[indice]
									resultmsms[0][5]=nombN[indice]
									resultmsms[0][6]=nombO[indice]
									resultmsms[0][7]=nombNa[indice]
									resultmsms[0][7]=W_sigma[1]
									resultmsms[0][8]=W_sigma[2]
									indiceresult+=1
								endif
							endif
						endif
					endif
				endif
				indice+=1
			while (massesimu(indice)>massemin)
		endif
	endif
	if (k15!=0)
//		if (abs(K15-resultmsms[0][1]))>0.1)
			if ((withD==0)&&(withNa==0))
				insertpoints indiceresult,1, resultmsms
				resultmsms[indiceresult][0]=bruitmoyen
				resultmsms[indiceresult][1]=K15
				resultmsms[indiceresult][2]=nombC[0]
				resultmsms[indiceresult][3]=nombH[0]
				resultmsms[indiceresult][4]=nombN[0]
				resultmsms[indiceresult][5]=nombO[0]
				resultmsms[indiceresult][6]=0
				resultmsms[indiceresult][7]=0
				indiceresult+=1
			else
				if (withNa==0)
					insertpoints indiceresult,1, resultmsms
					resultmsms[indiceresult][0]=bruitmoyen
					resultmsms[indiceresult][1]=K15
					resultmsms[indiceresult][2]=nombC[0]
					resultmsms[indiceresult][3]=nombH[0]
					resultmsms[indiceresult][4]=nombD[0]
					resultmsms[indiceresult][5]=nombN[0]
					resultmsms[indiceresult][6]=nombO[0]
					resultmsms[indiceresult][7]=0
					resultmsms[indiceresult][8]=0
					indiceresult+=1
				endif	
				if (withD==0)
					insertpoints indiceresult,1, resultmsms
					resultmsms[indiceresult][0]=bruitmoyen
					resultmsms[indiceresult][1]=K15
					resultmsms[indiceresult][2]=nombC[0]
					resultmsms[indiceresult][3]=nombH[0]
					resultmsms[indiceresult][4]=nombN[0]
					resultmsms[indiceresult][5]=nombO[0]
					resultmsms[indiceresult][6]=nombNa[0]
					resultmsms[indiceresult][7]=0
					resultmsms[indiceresult][8]=0
					indiceresult+=1
				endif
				if ((withD!=0)&&(withNa!=0))
					insertpoints indiceresult,1, resultmsms
					resultmsms[indiceresult][0]=bruitmoyen
					resultmsms[indiceresult][1]=K15
					resultmsms[indiceresult][2]=nombC[0]
					resultmsms[indiceresult][3]=nombH[0]
					resultmsms[indiceresult][4]=nombD[0]
					resultmsms[indiceresult][5]=nombN[0]
					resultmsms[indiceresult][6]=nombO[0]
					resultmsms[indiceresult][7]=nombNa[0]
					resultmsms[indiceresult][8]=0
					resultmsms[indiceresult][9]=0
					indiceresult+=1			
				endif	
			endif
//		endif
	endif
	deletepoints/M=0 indiceresult, 2, resultmsms
	make/O/N=(indiceresult), resultI, resultM
	resultI=resultmsms[x][0]
	resultM=resultmsms[x][1]
	dowindow/F evaluateur
	append resultI vs resultM
	ModifyGraph mode(resultI)=1,lsize(resultI)=3,rgb(resultI)=(65280,0,52224)
	indice=0
	string formule
	if ((withD==0)&&(withNa==0))
		formule="aaCHNO"
	else
		if (withNa==0)
			formule="aaCHSNO"
		endif
		if (withD==0)
			formule="aaCHNOn"
		endif
		if ((withD!=0)&&(withNa!=0))
			formule="aaCHSNOn"
		endif
	endif
	variable indicetag=2
	do
		letag+="\\JC"
		nomtag="text"+num2str(indice)
		do
			if (resultmsms[indice][indicetag]!=0)		
				letag+=formule[indicetag]
				if (resultmsms[indice][indicetag]>1)
					letag+="\\B"+num2str(resultmsms[indice][indicetag])+"\\M"
				endif
			endif
			indicetag+=1
		while (indicetag<=(5+withD+withNa))
		if (K9==-1)
			letag+="\\S-\\M\r"
		else
			letag+="\\S+\\M\r"
		endif
		letag+="\\K(65280,0,0)"
		indicetag=2
		do
			if ((resultmsms[indiceresult][indicetag]-resultmsms[indice][indicetag])!=0)		
				letag+=formule[indicetag]
				if ((resultmsms[indiceresult][indicetag]-resultmsms[indice][indicetag])>1)
					letag+="\\B"+num2str((resultmsms[indiceresult][indicetag]-resultmsms[indice][indicetag]))+"\\M"
				endif
			endif
			indicetag+=1
		while (indicetag<=(5+withD+withNa))
		Tag/C/N=$nomtag/F=0/X=-10/Y=0 /L=1/TL=0 resultI, indice ,letag
	indice+=1				
	letag=""
	indicetag=2
	while (indice<indiceresult)				
	zoom ()
end

Function zoomgauche(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	variable massebas
	variable massehaut
	nvar pointchek	
	variable point
	string tagname
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z lamassemsms
			wave/Z lintensite
			wave/Z resultI
			wave/Z resultM
			dowindow/F evaluateur
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65535,65535,65535) 
			point=pointchek-1
			pointchek=point
			if (pointchek>=0)
				massebas=(floor(resultM[pointchek]))* 0.9996
				massehaut=(floor(resultM[pointchek]))* 1.0016
				setaxis zoomaxe massebas, massehaut
			else 
				pointchek+=1
			endif
			execute "zoom()"
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65280,65280,0) 
			dowindow/F panel10
			break
		case 3: // mouse out
			wave/Z resultM
			dowindow/F evaluateur
			massebas=((resultM[pointchek]))* 0.99995
			massehaut=((resultM[pointchek]))* 1.00005
			setaxis zoomaxe massebas, massehaut
			break
	endswitch
	return 0
End

Function checkmsms(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z lamassemsms
			wave/Z lintensite
			wave/Z resultI
			wave/Z resultM
			variable massebas
			variable massehaut
			variable/G pointchek
			string result
			string tagname
			pointchek=numpnts(resultM)-1
			dowindow/F evaluateur
			result= axisinfo("","zoomaxe")
			DelayUpdate
			if (strlen(result)==0)
				AppendToGraph/B=zoomaxe lintensite vs lamassemsms
				AppendToGraph/B=zoomaxe resultI vs resultM
				ModifyGraph axisEnab(bottom)={0,0.65},axisEnab(zoomaxe)={0.66,1}
				ModifyGraph freePos(zoomaxe)=0
				ModifyGraph mode(resultI#1)=1,lsize(resultI#1)=1,rgb(resultI#1)=(65280,0,52224)
				AppendToGraph/L=simu/B=zoomaxe fakeint vs massesimu
				ModifyGraph mode(fakeint#1)=1,rgb(fakeint#1)=(0,0,0),rgb(lintensite#1)= (39168,26112,0)
				tagname = "text"+num2str(pointchek)
				Tag/C/N=$tagname/B=(65280,65280,0) 
			endif
			massebas=(floor(resultM[pointchek]))* 0.9996
			massehaut=(floor(resultM[pointchek]))* 1.0016
			setaxis zoomaxe massebas, massehaut
			dowindow/F panel10
			break
	endswitch
	return 0
End

Function zoomdroit(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z lamassemsms
			wave/Z lintensite
			wave/Z resultI
			wave/Z resultM
			variable point
			variable massebas
			variable massehaut
			nvar pointchek	
			string tagname
			dowindow/F evaluateur
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65535,65535,65535) 
			point=pointchek+1
			pointchek=point
			if (pointchek<(numpnts(resultI)))
				massebas=(floor(resultM[pointchek]))* 0.9996
				massehaut=(floor(resultM[pointchek]))* 1.0016
				setaxis zoomaxe massebas, massehaut
			else 
				pointchek-=1
			endif
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65280,65280,0) 
			execute "zoom()"
			dowindow/F panel10
			break
		case 3: // mouse out
			dowindow/F evaluateur
			massebas=((resultM[pointchek]))* 0.99995
			massehaut=((resultM[pointchek]))* 1.00005
			setaxis zoomaxe massebas, massehaut
			break
	endswitch
	return 0
End

Function killresult(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z lamassemsms
			wave/Z lintensite
			wave/Z resultI
			wave/Z resultM
			variable point
			nvar pointchek
			string tagname
			getaxis /W= evaluateur/Q left
			resultI[pointchek]=V_min
			deletepoints/M=0 pointchek, 1, resultmsms
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65535,65535,65535) 
			pointchek-=1
			variable massebas
			variable massehaut
			massebas=(floor(resultM[pointchek]))* 0.9996
			massehaut=(floor(resultM[pointchek]))* 1.0016
			tagname = "text"+num2str(pointchek)
			Tag/C/N=$tagname/B=(65280,65280,0) 
			setaxis zoomaxe massebas, massehaut
			execute "zoom()"
			dowindow/F panel10
			break
	endswitch
	return 0
End

Function savemsms(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			wave/Z resultmsms
			string nomsauvegarde
			svar chiffre
			nomsauvegarde="MSMSresult"+chiffre+".txt"
			Save/O/J/W/P=lasauvegarde resultmsms as nomsauvegarde
			GetSelection notebook,LoadedData,2
			notebook  LoadedData ruler=title, text= s_selection
			dowindow/F panel10
			break
	endswitch
	return 0
End
//==========MSMS ROLAND ENDS========
Function Wesslau(w,M) : FitFunc
	Wave w
	Variable M

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(m) = K* 1/M * 1/(sqrt(2*pi*N*p*(1-p))) * exp( -0.5* (ln(M/M0)/ln(a)-N*p)^2/ (N*p*(1-p)) )
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ M
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = N
	//CurveFitDialog/ w[1] = p
	//CurveFitDialog/ w[2] = M0
	//CurveFitDialog/ w[3] = a
	//CurveFitDialog/ w[4] = K

	return w[4]* 1/M * 1/(sqrt(2*pi*w[0]*w[1]*(1-w[1]))) * exp( -0.5* (ln(M/w[2])/ln(w[3])-w[0]*w[1])^2/ (w[0]*w[1]*(1-w[1])) )
End



// ROLAND  modifications fin janvier 2016

////////////////////////////////////////////////////////////////////////////////////////////////////////
// selectionFonctionNbO()
//	Fonction qui fait la selection en O. On va générer les waves de
//	selection en fonction du paramètre de slider. Si le choix est -1 la
//	selection va prendre l'ensemble des données.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function selectionFonctionNbO()
	NVar compareAttribNbOGraphs
	Wave element_O
	Wave MassInconnue
	
	FindLevel /P /Q element_O compareAttribNbOGraphs
	Variable debutSelection = V_levelX
	FindLevel /P /Q element_O compareAttribNbOGraphs+1
	Variable finSelection
	nvar charge
	If(V_Flag == 0)
		// Il y a le niveau n+1
		finSelection = V_LevelX - 1
	Else
		// Il n'y a pas le level n+1, on est intéréssé par le plus grand nombre d'O
		finSelection = numpnts(element_O) -1
	EndIf
	
	Duplicate /O /R=[debutSelection,finSelection] MassInconnue selection_masse
	
	// On creer les waves de selection
	If(waveexists(element_H)!=0)
		Duplicate /O /R=[debutSelection,finSelection] element_H selection_element_H
	Else
		// on crée la selection azote que l'on rempli de 0
		Make/O/N=(finSelection-debutSelection) selection_element_H = 0
	EndIf
	
	If(waveexists(element_N)!=0)
		Duplicate /O /R=[debutSelection,finSelection] element_N selection_element_N
	Else
		// on crée la selection azote que l'on rempli de 0
		Make/O/N=(finSelection-debutSelection) selection_element_N = 0
	EndIf
	
	If(waveexists(element_C)!=0)
		Duplicate /O /R=[debutSelection,finSelection] element_C selection_element_C
	Else
		// on crée la selection azote que l'on rempli de 0
		Make/O/N=(finSelection-debutSelection) selection_element_C = 0
	EndIf
	
	If(waveexists(element_O)!=0)
		Duplicate /O /R=[debutSelection,finSelection] element_O selection_element_O
	Else
		// on crée la selection azote que l'on rempli de 0
		Make/O/N=(finSelection-debutSelection) selection_element_O = 0
	EndIf
	
	// Calcul des fonctions utilisées dans graphs
	Duplicate /O selection_element_H selection_HsurC
	selection_HsurC = ((selection_element_H-charge) - selection_element_N) / selection_element_C
	Duplicate /O selection_element_N selection_NsurC
	selection_NsurC = selection_element_N / selection_element_C
	Duplicate /O selection_element_C selection_DBE
	selection_DBE = selection_DBE - ((selection_element_H-charge)/2) + (selection_element_N/2) +1
	Duplicate /O selection_element_C selection_excesC
	selection_excesC = selection_excesC - (((selection_element_H-charge)+selection_element_N)/2)
	Duplicate /O selection_element_H selection_excesMethylene
	selection_excesMethylene = ((selection_excesMethylene-charge) - (3*selection_element_N)) / 2
	Duplicate /O selection_excesC selection_composite
	selection_composite = selection_composite + (0.5*selection_excesMethylene/25)
	
	If(Wintype("GraphsAttributionComparaison#affichageLegendColorDBE")!=0)
		//Mise a jour de l'echelle de couleur
		ColorScale /C /N=colorScaleDBE /W=GraphsAttributionComparaison#affichageLegendColorDBE ctab={wavemin(selection_DBE),wavemax(selection_DBE),CyanMagenta}
	endIf
End // selectionFonctionNbO()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// interfaceGraphsCompAttrib()
//	Création de l'interface affichant Graphs, avec selection en nombre
//	d'oxygène.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function interfaceGraphsCompAttrib()
	
	Wave /T compareEnCoursElementsWaves
	Wave MassInconnue
	
	// Variables permettent de choisir combien on veut de graph par ligne et par colonne... surement inutile mais au cas où l'option est présente.
	Variable nbGraphsCollonesCompareAttrib = 4
	Variable	nbGraphsLignesCompareAttrib = 2
	
	If(Wintype("GraphsAttributionComparaison")!=0)
		KillWindow GraphsAttributionComparaison
	endIf
	
	// On fait un panneau de taille liée à l'écran principal
	NewPanel /N=GraphsAttributionComparaison /W=(10,60,extractionLargueurEcran()-20,extractionHauteurEcran()-120)
	SetWindow GraphsAttributionComparaison ,hook(hookAttributionCompareGraphs)=hookAttributionCompareGraphs,hookevents=1
	
	// Affichage du nom de la molécule actuellment sélectionée
	TitleBox selectionMoleculeGraphs win=GraphsAttributionComparaison
	TitleBox selectionMoleculeGraphs pos={760,2},size={300,31}
	TitleBox selectionMoleculeGraphs fixedSize=1
	TitleBox selectionMoleculeGraphs fStyle=1
	TitleBox selectionMoleculeGraphs frame=1
	TitleBox selectionMoleculeGraphs anchor=LC
	TitleBox selectionMoleculeGraphs title="Molécule : NONE"
	
	// Mise en place si necessaire du slider de tri en Oxygène
	// variable globale va contenir le nombre d'oxygène désiré, si = -1 -> pas de selection en oxygène
	Variable /G compareAttribNbOGraphs = 0
	Variable minOAttrib = 0
	Variable maxOAttrib = 0
	If(waveexists(element_O)!=0)
		minOAttrib = wavemin(element_O)
		compareAttribNbOGraphs = minOAttrib-1
		maxOAttrib = wavemax(element_O)
		
		TitleBox titreSlider win=GraphsAttributionComparaison,pos={5,5},size={30,20},title="Nb O",fsize=15,fstyle=1,frame=0
		Slider reglageOCompareAttrib win=GraphsAttributionComparaison,pos={0,35},size={5,250},limits={minOAttrib-1,maxOAttrib,1},side=2,vert=1,ticks=(maxOAttrib-minOAttrib)+1,fsize=1,variable=compareAttribNbOGraphs,proc=selectionNbOGraphsCompareAttrib
		TitleBox legendSlider win=GraphsAttributionComparaison,pos={2,265},size={15,25},title="All",fsize=12,fstyle=1,frame=0,labelBack=(60000,60000,60000),fixedSize=1
		
		// fait le tri en O
		Variable j
		Variable nombreElements = numpnts(compareEnCoursElementsWaves)
		Duplicate /O element_O elementTriOLD
		For(j=0 ; j<nombreElements ; j+=1)
			Sort elementTriOLD,$compareEnCoursElementsWaves[j]
		endFor
		Sort elementTriOLD,Delta,Proximite,Intensitepic,logintensitePic,MassInconnue,listbagCompareEnCours
		KillWaves /Z elementTriOLD
		
		// On a fait un tri en Oxygène, on bloque le menu de tri dans l'interface principale de Compare Attrib
		// pour éviter que les données se retrouvent toutes remélangées.
		PopupMenu choixTri win=AttributionComparaison,popmatch="O"
		PopupMenu choixTri win=AttributionComparaison,disable=2
	
		// On fait la selection en Oxygènhe avec  nbO = 0
		selectionFonctionNbO()	
	Else
		
		// Pas d'oxygene dans l'attribution
		// on va créer les waves de selection necessaires pour l'affichage des graphs basés sur des selections.
		// on place la wave entière dans la wave de selection. Si l'element n'existe pas alors on place une wave vide pour que celle-ci existe.
		Duplicate /O MassInconnue selection_masse
		If(waveexists(element_H)!=0)
			Duplicate /O element_H selection_element_H
		Else
			Make/O/N=(numpnts(element_H)) selection_element_H = 0
		endIf
		If(waveexists(element_N)!=0)
			Duplicate /O element_N selection_element_N
		Else
			Make/O/N=(numpnts(element_N)) selection_element_N = 0
		endIf
		If(waveexists(element_C)!=0)
			Duplicate /O element_C selection_element_C
		Else
			Make/O/N=(numpnts(element_C)) selection_element_C = 0
		endIf
		If(waveexists(element_O)!=0)
			Duplicate /O element_O selection_element_O
		Else
			Make/O/N=(numpnts(element_O)) selection_element_O = 0
		endIf
	endIf
	
	// Dans tout les cas, présence d'oxygène ou non, on calcule les waves non liées au sélection mais au jeu de donnée entier.
	Wave element_H
	Wave element_C
	Wave element_N
	Wave element_O
	Duplicate /O element_C Zubarev
	Wave Zubarev
	Zubarev = Zubarev - ((element_N + element_H -1) /2)
	
	Duplicate /O element_O VanKrev_OsurC
	VanKrev_OsurC = VanKrev_OsurC / element_C
	Wave VanKrev_OsurC
	
	Duplicate /O element_H VanKrev_HsurC
	VanKrev_HsurC = ((VanKrev_HsurC-1) - element_N) / element_C
	Wave VanKrev_HsurC
	
	Duplicate /O element_N VanKrev_NsurC
	VanKrev_NsurC = VanKrev_NsurC / element_C
	Wave VanKrev_NsurC
	
	Duplicate /O MassInconnue VanKrevMass
	
	// Contruire les graphiques
	// declare waves de selection
	Wave selection_HsurC
	Wave selection_NsurC
	Wave selection_DBE
	Wave selection_excesC // apellé aussi le Zubarev
	Wave selection_excesMethylene
	Wave selection_composite

	Wave selection_masse
	
	// Recup de la zone disponible pour graphs
	Variable largueurDispo = extractionLargueurEcran()-85
	Variable largueurGraph = largueurDispo / nbGraphsCollonesCompareAttrib
	Variable hauteurDispo = extractionHauteurEcran()-220
	Variable hauteurGraph = hauteurDispo / nbGraphsLignesCompareAttrib
	
	// mise en place des graphs
	Display /HOST=GraphsAttributionComparaison /N=graph11 /W=(45+0*largueurGraph , 35+0*hauteurGraph , 45+1*largueurGraph , 35+1*hauteurGraph) selection_HsurC vs selection_masse
	ModifyGraph /W=GraphsAttributionComparaison#graph11 gbRGB=(30000,20000,10000),mode[0]=3,marker[0]=13,msize[0]=4,mrkThick[0]=1.5,zColor[0]={selection_DBE,*,*,CyanMagenta}
	Label left "(H-N)/C"
	Label bottom "m/z"
	
	Display /HOST=GraphsAttributionComparaison /N=graph21 /W=(45+0*largueurGraph , 35+1*hauteurGraph , 45+1*largueurGraph , 35+2*hauteurGraph) selection_NsurC vs selection_masse
	ModifyGraph /W=GraphsAttributionComparaison#graph21 gbRGB=(30000,20000,10000),mode[0]=3,marker[0]=13,msize[0]=4,mrkThick[0]=1.5,zColor[0]={selection_DBE,*,*,CyanMagenta}
	Label left "N/C"
	Label bottom "m/z"
	
	Display /HOST=GraphsAttributionComparaison /N=graph12 /W=(45+1*largueurGraph , 35+0*hauteurGraph , 45+2*largueurGraph , 35+1*hauteurGraph) selection_DBE vs selection_masse
	ModifyGraph /W=GraphsAttributionComparaison#graph12 gbRGB=(30000,20000,10000),mode[0]=3, marker=41,msize[0]=4,mrkThick[0]=1.5,zColor[0]={selection_DBE,*,*,CyanMagenta}
	Label left "DBE=C-H/2+N/2+1"
	Label bottom "m/z"
	
	Display /HOST=GraphsAttributionComparaison /N=graph22 /W=(45+1*largueurGraph , 35+1*hauteurGraph , 45+2*largueurGraph , 35+2*hauteurGraph) selection_composite vs selection_masse
	ModifyGraph /W=GraphsAttributionComparaison#graph22 gbRGB=(30000,20000,10000),mode[0]=3,marker[0]=8,mrkThick[0]=2,zColor[0]={selection_DBE,*,*,CyanMagenta}, zero(left)=1
	Label left "Composite"
	Label bottom "m/z"
	
	Display /HOST=GraphsAttributionComparaison /N=graph13 /W=(45+2*largueurGraph , 35+0*hauteurGraph , 45+3*largueurGraph , 35+1*hauteurGraph) Zubarev vs VanKrevMass
	ModifyGraph /W=GraphsAttributionComparaison#graph13 gbRGB=(40000,36000,40000),mode[0]=3,marker[0]=43,mrkThick[0]=1.5,zero(left)=1
	Label left "Zubarev = C - (N + H)/2"
	Label bottom "m/z"
	If(wavemax(element_O)>0)
		ModifyGraph /W=GraphsAttributionComparaison#graph13 zColor[0]={element_O,*,*,BlueHot}
	endIf
	
	Display /HOST=GraphsAttributionComparaison /N=graph23 /W=(45+2*largueurGraph , 35+1*hauteurGraph , 45+3*largueurGraph , 35+2*hauteurGraph) selection_excesMethylene vs selection_excesC
	ModifyGraph  /W=GraphsAttributionComparaison#graph23 gbRGB=(40000,36000,40000),mode[0]=3,marker[0]=8,mrkThick[0]=2,zmrkSize[0]={selection_masse,*,*,15,1},zColor[0]={selection_masse,0,*,YellowHot},zero=1
	Label left "Methylene excess = (H-3N)/2"
	Label bottom "Carbon excess = Zubarev = C-(N+H)/2"
	
	Display /HOST=GraphsAttributionComparaison /N=graph14 /W=(45+3*largueurGraph , 35+0*hauteurGraph , 45+4*largueurGraph , 35+1*hauteurGraph) VanKrev_HsurC vs VanKrev_OsurC
	ModifyGraph /W=GraphsAttributionComparaison#graph14 gbRGB=(40000,36000,40000), mode[0]=3,marker[0]=8,mrkThick[0]=2,zmrkSize[0]={VanKrevMass,*,*,8,0.5}
	Label left "(H-N)/C"
	Label bottom "O/C"
	If (wavemax(element_O)>0)
		ModifyGraph /W=GraphsAttributionComparaison#graph14 zColor[0]={element_O,*,*,BlueHot}
	endIf
	// ajout lignes Vankrev
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,2,0.204129604559247,0.56232429042305
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.251804610887229,0.484848499298096
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.313916414906079,0.422368022584424
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.358905505384598,0.554826633217408
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.427732099027107,0.714776653604409
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.546919614847062,0.907216521882519
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.580493562965359,1.41705721186608
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.582507999852457,2.00187447390606
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.5,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.559677715132015,1.67197755685787
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.572435815416968,1.15213999060011
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.438475762424962,2.18681668497852
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.371999345150734,2.22930340914382
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.276985071975952,2.20930965659545
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.217894923287749,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.182306538282354,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.156790337712448,2.18181824684143
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.138996145209751,2.21930653286963
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,0,0.121201952707053,2.17931902777288
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.185328193613001,0.512339909052111
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.158469035118363,0.412371146310236
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.137653187285019,0.339893793322376
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.122880650112968,0.274914097540157
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1.99937525483751,0.109115331384466,0.247422687786142
	
	
	Display /HOST=GraphsAttributionComparaison /N=graph24 /W=(45+3*largueurGraph , 35+1*hauteurGraph , 45+4*largueurGraph , 35+2*hauteurGraph) VanKrev_HsurC vs VanKrev_NsurC
	ModifyGraph /W=GraphsAttributionComparaison#graph24 gbRGB=(40000,36000,40000),mode[0]=3,marker[0]=8,mrkThick[0]=2,zmrkSize[0]={VanKrevMass,*,*,8,0.5},zColor[0]={VanKrevMass,0,*,YellowHot}
	Label left "(H-N)/C"
	Label bottom "N/C"
	//SetAxis left 0,2.6
	SetAxis bottom 0,2.6
	// ajout lignes Vankrev
	SetDrawLayer UserFront
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1,2
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 0,1,2,1
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.22820211515864,0.923253150057274
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.48025851938895,0.956013745704467
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.441480611045828,1.67079037800687
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.68954171562867,0.920274914089346
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.24500587544066,1.45635738831615
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,1.18848413631022,1.49805269186712
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.742420681551117,1.54868270332188
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.543830787309048,1.82565864833906
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.661457109283196,1.71546391752577
	SetDrawEnv xcoord= bottom,ycoord= left
	DrawLine 1,0,0.378848413631022,2.08178694158076
	
	
	// On ajoutte les echelles de couleurs qui seront mise à jour pour suivre l'évolution du jeu de donnée.
	Display /HOST=GraphsAttributionComparaison /N=affichageLegendColorDBE /W=(0,360,45,660)
	ColorScale /W=GraphsAttributionComparaison#affichageLegendColorDBE /N=colorScaleDBE heightPct=84,vert=1,tickLen=-5,ctab={wavemin(selection_DBE),wavemax(selection_DBE),CyanMagenta}
	TitleBox titreDBE win=GraphsAttributionComparaison,pos={2,340},size={20,20},title="DBE",fsize=15,fstyle=1,frame=0
	
	Display /HOST=GraphsAttributionComparaison /N=affichageLegendColorMass /W=(450,0,750,40)
	ColorScale /W=GraphsAttributionComparaison#affichageLegendColorMass /A=LC /N=colorScaleMass widthPct=90,vert=0,tickLen=-5,ctab={wavemin(VanKrevMass),wavemax(VanKrevMass),YellowHot}
	TitleBox titreMass win=GraphsAttributionComparaison,pos={405,5},size={20,20},title="Masse",fsize=15,fstyle=1,frame=0
	
	Display /HOST=GraphsAttributionComparaison /N=affichageLegendColorElementO /W=(45,0,345,40)
	ColorScale /W=GraphsAttributionComparaison#affichageLegendColorElementO /N=colorScaleElementO widthPct=90,tickLen=-5,vert=0,ctab={wavemin(element_O),wavemax(element_O),BlueHot}
End // interfaceGraphsCompAttrib()

////////////////////////////////////////////////////////////////////////////////////////////////////////
// supprimeMoleculeCompareAttrib(nomBouton)
//	Fonction supprimant la molécule actuellement sélectionée
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function supprimeMoleculeCompareAttrib(nomBouton) : ButtonControl
	String nomBouton
	
	NVar pointeurAttribCompareSelec
	
	Wave/T listbagCompareEnCours
	Wave/T compareEnCoursElementsWaves
	wave MassInconnue,intensitePic, logintensitePic
	wave gicm=massegiclee
	wave gicint=intensitegiclee
	insertpoints 0,1, gicm, gicint
	gicm[0]=MassInconnue[pointeurAttribCompareSelec]
	gicint[0]=intensitePic[pointeurAttribCompareSelec]
	DeletePoints pointeurAttribCompareSelec,1,MassInconnue, Delta, intensitePic,logintensitePic, Proximite,listbagCompareEnCours
     	goods2data("bel");gicle2data("gic");  

	Variable j
	For(j=0 ; j < numpnts(compareEnCoursElementsWaves) ; j+=1)
		DeletePoints pointeurAttribCompareSelec,1,$compareEnCoursElementsWaves[j]
	endFor
	
End // supprimeMoleculeCompareAttrib(nomBouton)

////////////////////////////////////////////////////////////////////////////////////////////////////////
//  interfaceCompareAttribution()
//	Créer l'interface de comparaison des attributions.
////////////////////////////////////////////////////////////////////////////////////////////////////////
Function interfaceCompareAttribution()
	Wave /T compareEnCoursElementsWaves
	Wave /T molbagCompareEnCoursElements
	Variable nombreElements = numpnts(compareEnCoursElementsWaves)
	
	Wave Delta
	Wave Proximite
	Wave MassInconnue
	Wave Intensitepic
	duplicate/O intensitepic logintensitepic
	logintensitepic= log(intensitepic)
	
	NewPanel /K=1/N=AttributionComparaison /W=(100,80,1500,780)
	SetWindow AttributionComparaison ,hook(hookAttributionCompare)=hookAttributionCompare,hookevents=1
	
	//1. Graph de d'écart PPM
	Display /HOST=AttributionComparaison /N=DeltaMass /W=(5,60,1395,155) Delta vs MassInconnue
	ModifyGraph /W=AttributionComparaison#DeltaMass mode[0] = 3,grid(left)=1, zero(left)=1, axThick(bottom)=0, margin(bottom)=8, margin(top)=8,margin(left)=70, margin(right)=10, rgb[0]=(20000,20000,50000),lblMargin(left)=5
	Label/W=AttributionComparaison#DeltaMass left "\\Z12Mass Bias\r(PPM)"
	
	// 2. Graph de proximité
	Display /HOST=AttributionComparaison /N=Proximite /W=(5,155,1395,235) Proximite vs MassInconnue
	ModifyGraph /W=AttributionComparaison#Proximite mode[0] = 1,grid(left)=1, log(left)=1, axThick(bottom)=0, margin(bottom)=8, margin(top)=8, margin(left)=70, margin(right)=10,rgb[0]=(0,0,0),lblMargin(left)=5
	Label/W=AttributionComparaison#Proximite left "\\Z12Loss function"
	
	// 3 Graphs du nombre d'atome vs masse
	Display /HOST=AttributionComparaison /N=nombreAtome /W=(5,235,1395,695) $compareEnCoursElementsWaves[0] vs MassInconnue
	ModifyGraph /W=AttributionComparaison#nombreAtome margin(top)=8, margin(left)=70, margin(right)=10,lblMargin(left)=5
	ModifyGraph /W=AttributionComparaison#nombreAtome mode[0]=3, marker[0]=41,zmrkSize[0]={logintensitepic,*,*,1,10},rgb[0]=(0,0,0)
	Label/W=AttributionComparaison#nombreAtome left "\\Z14Nb atomes"
	Label/W=AttributionComparaison#nombreAtome bottom "\\Z16Masse/charge"
	// ajoute les autres éléments de façon réccursive en créant la légende
	String textLengend = "\s(#0) "+ molbagCompareEnCoursElements[0] + " "
	String nomCheckBoxCreation = "affichageElement0"
	CheckBox $nomCheckBoxCreation win=AttributionComparaison,title=molbagCompareEnCoursElements[0],size={40,15},pos={1055,5},value=1,proc=affichageElementChanged
	Wave rougeCompareAttrib
	Wave vertCompareAttrib
	Wave bleuCompareAttrib
	Variable j
	For(j=1 ; j<numpnts(compareEnCoursElementsWaves) ; j+=1)
		AppendToGraph /W=AttributionComparaison#nombreAtome $compareEnCoursElementsWaves[j] vs MassInconnue
		ModifyGraph /W=AttributionComparaison#nombreAtome mode[j]=3, marker[j]=41,mrkThick[j]=1.5,rgb[j]=(rougeCompareAttrib[j],vertCompareAttrib[j],bleuCompareAttrib[j])
		nomCheckBoxCreation = "affichageElement" + num2str(j)
		CheckBox $nomCheckBoxCreation win=AttributionComparaison, title=molbagCompareEnCoursElements[j],size={35,15},pos={1055 + (j-mod(j,3))/3 *40 , 5+ (mod(j,3)) * 18},value=1,proc=affichageElementChanged
		textLengend = textLengend + "\r\s(#" + num2str(j) + ") "+ molbagCompareEnCoursElements[j] + " "
	endFor
	Legend /W=AttributionComparaison#nombreAtome /N=legendNombreAtome /A=LT textLengend
	
	// 5. Mise en place de l'interface
	// Choix molbag et sauvegarde
	string camembert
	PopupMenu choixMolBag win=AttributionComparaison,title="MolBag chargé ", value=replacestring("Molbag_",wavelist("Molbag_*",";",""),""),pos={5,5},size={245,25},bodyWidth=165,proc=chargerMolbagAttribCompare
	PopupMenu choixMolBag popmatch="Molbag_Current"
	
	TitleBox sep1 win=AttributionComparaison,pos={255,0},size={5,60},fixedSize=1,frame=4,title=" "
	
	// Suppression de la molécule sélectionée
	Button supprimeMolecule win=AttributionComparaison,title="Supprime Molécule", pos={265,5},size={175,25},fColor=(54000,3600,2300),fStyle=1,proc=supprimeMoleculeCompareAttrib
	
	// Menu de tri, et placement sur le tri en Procimité qui est celui par défaut de l'attribution.
	Camembert =genereListeFromTextWave(molbagCompareEnCoursElements)
	PopupMenu choixTri win=AttributionComparaison,title="Tri en fonction de ", value=#camembert,pos={255,35},size={180,25},bodyWidth=70,proc=trierAttribCompare
	PopupMenu choixTri popmatch="Proximité"
	
	// Boutons déplacement dans les données
	Button precAttribCompare win=AttributionComparaison,title="Prec", pos={450,5},size={35,50},fColor=(65000,50000,3600),fStyle=1,proc=prevCompareAttribBut
	Button suivAttribCompare win=AttributionComparaison,title="Suiv", pos={485,5},size={35,50},fColor=(65000,50000,3600),fStyle=1,proc=nextCompareAttribBut
	
	// Affichage nom de la molécule sélectionée
	TitleBox selectionMolecule win=AttributionComparaison
	TitleBox selectionMolecule pos={530,10},size={300,40}
	TitleBox selectionMolecule fixedSize=1
	TitleBox selectionMolecule fStyle=1
	TitleBox selectionMolecule frame=1
	TitleBox selectionMolecule anchor=LC
	TitleBox selectionMolecule title="Molécule : NONE"
	// Affichage infos molécule sélectionée
	TitleBox affichageInfosMolecule1 win=AttributionComparaison
	TitleBox affichageInfosMolecule1 pos={840,5},size={150,50}
	TitleBox affichageInfosMolecule1 fixedSize=1
	TitleBox affichageInfosMolecule1 fStyle=1
	TitleBox affichageInfosMolecule1 frame=0
	TitleBox affichageInfosMolecule1 anchor=LC
	TitleBox affichageInfosMolecule1 title="Masse : NONE\rIntensite : NONE\rProximite : NONE\rDelta : NONE"
	
	TitleBox sep2 win=AttributionComparaison,pos={975,0},size={5,60},fixedSize=1,frame=4,title=" "
	
	Button timeForGraph win=AttributionComparaison,title="Graphs", pos={990,5},size={50,50},fColor=(8700,45000,20000),fStyle=1,proc=graphsCompareAttrib
	
	TitleBox sep3 win=AttributionComparaison,pos={1045,0},size={5,60},fixedSize=1,frame=4,title=" "

	// Place les curseurs sur les graphs
	changeSelectionAttribCompare()
End // interfaceCompareAttribution()

Function ListBoxProcDataMan(lba) : ListBoxControl
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
	string lax, lay, laz
	NVar spectrummode

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
		if(row==-1)
			resortDATAbycol(selwave,listwave,col)
		endif
			break
		case 3: // double click
			break
		case 4: // cell selection
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			if(col==3)
				lay=listwave[row][1]+"_1stripped"
				lax=listwave[row][1]+"_0stripped"
				switch(selwave[row][col])	// numeric switch
					case 48:		// execute if case matches expression
						AppendToGraph/W=agregateur $lay vs $lax
						ModifyGraph/W=agregateur rgb($lay)=(60000*(ln(stringcrc(0,listwave[row][1]))-18)/5,str2num(listwave[row][2]),100*row)
						ModifyGraph/W=agregateur mode($lay)=spectrummode
						setaxis/W=agregateur/A
						break						// exit from switch
					case 32:		// execute if case matches expression
						RemoveFromGraph/W=agregateur $lay
						break
				endswitch
			elseif(col==4)
				lay=listwave[row][1]+"_1stripped"
				lax=listwave[row][1]+"_0stripped"
				laz=listwave[row][1]+"_2stripped"
				switch(selwave[row][col])	// numeric switch
					case 48:		// execute if case matches expression
						wave transx=$lax
						duplicate/O $lay trans
						trans=defmass(transx)
						duplicate/O trans $laz
						killwaves trans
						AppendToGraph/W=dmvm $laz vs $lax
						ModifyGraph/W=dmvm rgb($laz)=(60000*(ln(stringcrc(0,listwave[row][1]))-18)/5,str2num(listwave[row][2]),100*row)
						ModifyGraph/W=dmvm mode($laz)=3,marker($laz)=19, msize($laz)=1.5
						//setaxis/W=dmvm/A
						break						// exit from switch
					case 32:		// execute if case matches expression
						RemoveFromGraph/W=dmvm $laz
						killwaves $laz
						break
				endswitch
			endif
			break
	endswitch

	return 0
End

Window filtres() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(436,651,1781,685) as "Filtres"
	Button Bouton1,pos={103,2},size={92,30},proc=bouton1,title="Union Find + radio"
	Button Bouton2,pos={355,2},size={65,30},proc=bouton2,title="anti-Ringing"
	Button Bouton3,pos={420,2},size={90,30},proc=ButtonProc_3,title="Hi-Res Calibration"
	Button Bouton3,fColor=(32768,65280,0)
	Button button0,pos={657,2},size={60,30},proc=ButtonProc_5,title="ROI 2 data"
	Button button0,fColor=(65280,43520,0)
	Button button1,pos={601,2},size={55,30},proc=ButtonProc_6,title="Invert ROI"
	Button button1,fColor=(65280,43520,0)
	Button buttonFatNoise,pos={197,2},size={60,16},proc=ButtonKillFatNoise,title="Fat Noise"
	Button buttonFatNoise,help={"Remove the less intense points so the new intensity distribution is strictly convex"}
	Button buttonFatNoise,fColor=(65280,43520,0)
	Button buttonFineNoise,pos={197,17},size={60,16},proc=ButtonKillFineNoise,title="Fine Cut"
	Button buttonFineNoise,help={"Deletes the less intense points so the most intense peaks becomes twice less frequent."}
	Button buttonFineNoise,fColor=(65280,43520,0)
	Button buttonNormalize,pos={259,2},size={32,30},proc=ButtonNorm,title="Norm"
	Button buttonNormalize,help={"Divides the current sample data by its maximum intensity"}
	Button buttonNormalize,fColor=(65280,43520,0)
	SetVariable setvar0,pos={3,1},size={100,16},proc=updateUF,title="Union-Find"
	SetVariable setvar0,limits={0,inf,0.1},value= _NUM:0.8
	SetVariable setvar1,pos={3,16},size={100,16},proc=UpdateRadio,title="Radio"
	SetVariable setvar1,limits={0,inf,1},value= _NUM:30
	Button button2,pos={294,2},size={60,30},proc=ButtonProc_12,title="Convolver"
	Button button2,fColor=(32768,65280,0)
	Button button3,pos={715,2},size={90,30},proc=ButtonProc_15,title="Crop Ticked Data"
	Button button3,fColor=(32768,65280,0)
	Button Bouton5,pos={511,2},size={90,30},proc=ButtonProc_16,title="find same peaks"
	Button Bouton5,help={"will put in ROI the peaks that appear both in current and chosen datasets within ginven uncertainty"}
	Button Bouton5,fColor=(32768,65280,0)
EndMacro

Function ButtonProc_16(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			findsame()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

// de Aurélien Fresneau mars 2016

function findsame()
	string tocompare
	Variable ppm =5
	Prompt tocompare, "choose wave to compare: ",  popup,waveList("*_0stripped",";","")	
	Prompt ppm, "define precision: "		
	DoPrompt "Enter info for comparison", tocompare, ppm
	if (V_Flag)
		return -1								// User canceled
	endif
	wave target = $tocompare
    	wave/Z wave0stripped, wave2stripped, wave1stripped
   	 variable k, j, n=numpnts(wave0stripped), m=numpnts(target)
    	make/O/N=0 roipnts
   	 for(k=0;k<n;k+=1)
     	  	for(j=0;j<m;j+=1)
            		if(abs(wave0stripped[k]-target[j])<(ppm*1E-6)*wave0stripped[k])
                	insertpoints 0,1, roipnts
                	roipnts[0]=k
            		endif
        	endfor
    	endfor
   	duplicate/O roipnts roi1 roi0 roi2
    	roi1=wave1stripped[roipnts[x]]
    	roi0=wave0stripped[roipnts[x]]
    	roi2=wave2stripped[roipnts[x]]
end

function showisotopes(numeroelement)
	variable numeroelement
	wave abu=mendev_abundances
	wave mas=mendev_masses, molecule
	wave/T nom=elements
	make/O/N=10 showmassexact, showintexact,defshowmassexact, logdeshow
	showintexact=abu[numeroelement][x]
	logdeshow=log(showintexact)
	showmassexact=mas[numeroelement][x]
	defshowmassexact=showmassexact-round(showmassexact)
	sort /R showintexact, logdeshow, showmassexact, defshowmassexact, showintexact
	findlevel/Q  showintexact, 0
	deletepoints v_levelx, 10, showintexact, showmassexact,defshowmassexact,logdeshow
	wavestats/Q defshowmassexact
	variable moyenne=V_avg
	killisowin()
	display/W=(35.25,41.75,438,357.5)/N=isotopes showintexact vs showmassexact
	appendtograph/L=def defshowmassexact vs showmassexact
	ModifyGraph mode(showintexact)=1,lsize(showintexact)=2
	ModifyGraph log(left)=1
	SetAxis/A=2 left 0.0001,*
	SetAxis def moyenne-0.01,moyenne+0.01
	TextBox/N=text0/F=0/A=MC/X=36/Y=41  "\\Z24"+nom[numeroelement]
	ModifyGraph lblPos(left)=78,lblPos(def)=75
	ModifyGraph axisEnab(left)={0,0.49}
	ModifyGraph axisEnab(def)={0.52,1}
	ModifyGraph mode(defshowmassexact)=3, marker(defshowmassexact)=19
	ModifyGraph zero(def)=1
	ModifyGraph freePos(def)=0
	ModifyGraph zmrkSize(defshowmassexact)={logdeshow,-4,0,1,10}
	ModifyGraph zColor(defshowmassexact)={logdeshow,-4,0,YellowHot}
	Label bottom "\\Z14m/z"
	Label left "\\Z14Intensity"
	Label def "\\Z14Mass defect"
end

function killisowin()
string lesisowin=winlist("isotopes*",";",""),lawin
variable k, n=itemsinlist(lesisowin)
	for(k=0;k<n;k+=1)	// Initialize variables;continue test
		lawin=stringfromlist(k,lesisowin)
		killwindow $lawin
	endfor												// Execute body code until continue test is FALSE
end

Function chargerMolbagFacteurAff(nomMenu,indexSelect,strSelect) : PopupMenuControl
	String nomMenu
	Variable indexSelect
	String strSelect
	
	// Va charger dans  les masse des facteurs a afficher contenu dans le molbag selectioné listeMasseFacteursRosePentes
	// Va charger la liste correspondate dans  listeFacteursRosePentes
	
	Wave listeMasseFacteursRosePentes
	Wave/T listeFacteursRosePentes
	// on va dans un premier temps extraire le listBag Correspondant
	String nomListeBag = "Listbag_" + strSelect
	String nomMolbag = "Molbag_"+ strSelect
	Duplicate /O /T $nomListeBag listeFacteursRosePentes
	// On va maintenant génere la wave des masses correspondantes a ces facteurs
	Duplicate /O /D $nomMolbag molBagExtractionMassesEnCours
	Variable nombreMolecule = DimSize(molBagExtractionMassesEnCours,1)
	Make /O /N=(nombreMolecule) masseExtractionMLolbag
	Variable o = 0
	For(o=0 ; o<nombreMolecule ; o+=1)
		// On duplique la Nieme colone contenant la compo de la molécule
		Duplicate /O /D /R=[0,113][o,o] molBagExtractionMassesEnCours moleculeCalculMasseEnCours
		// On la multiplie par les masse d'élements
		Wave mendev_masses
		moleculeCalculMasseEnCours=moleculeCalculMasseEnCours*mendev_masses
		// On recupere la somme de tout qui correspond à la masse de la molécule		
		masseExtractionMLolbag[o] = sum(moleculeCalculMasseEnCours)
	endFor
	Duplicate /O masseExtractionMLolbag listeMasseFacteursRosePentes
	
	If(Wintype("RoseDesPentes")!=0)
		miseAJourRosePentes()
	else
		miseajourHTmarks()
	endIf
End // chargerMolbagFacteurAff(nomMenu,indexSelect,strSelect)


function miseajourHTmarks()
	Wave listeMasseFacteursRosePentes
	wave/T listefacteursrosepentes
	variable j
	string letag
	nvar borneHoughSlopes2, borneHoughSlopes1, borneHoughIntercepts2, borneHoughIntercepts1, storageHoughResolution,intervalleHoughSlopes,intervalleHoughIntercepts
	Duplicate/O listeMasseFacteursRosePentes listeSlopesFacteursHT
	listeSlopesFacteursHT= defmass(listeMasseFacteursRosePentes)/listeMasseFacteursRosePentes
	Make/O/D/N=0 intermarks, slomarks
		For(j=0 ; j< numpnts(listeMasseFacteursRosePentes) ; j+=1)
			letag="text"+num2str(j)
			if(listeSlopesFacteursHT[j]<borneHoughSlopes2 && listeSlopesFacteursHT[j]>borneHoughSlopes1)
				insertpoints numpnts(intermarks),3, intermarks, slomarks
				intermarks[numpnts(intermarks)-3]=(borneHoughIntercepts2-borneHoughIntercepts1)*(storageHoughResolution/intervalleHoughIntercepts)
				intermarks[numpnts(intermarks)-2]=(borneHoughIntercepts1-borneHoughIntercepts1)*(storageHoughResolution/intervalleHoughIntercepts)
				slomarks[numpnts(intermarks)-3]=(listeSlopesFacteursHT[j]-borneHoughSlopes1)*(storageHoughResolution/intervalleHoughSlopes)
				slomarks[numpnts(intermarks)-2]=(listeSlopesFacteursHT[j]-borneHoughSlopes1)*(storageHoughResolution/intervalleHoughSlopes)
				slomarks[numpnts(intermarks)-1]=NaN
				intermarks[numpnts(intermarks)-1]=NaN
				Tag/W=etudeTransformationHough#imageResultatHough/C/N=$letag/L=1/TL={lineRGB=(65535,65535,65535)}/Y=-10 intermarks, numpnts(intermarks)-3, listefacteursrosepentes[j]
			else
				Tag/W=etudeTransformationHough#imageResultatHough/K/N=$letag
			endif
		endFor
end

function pix2inter(pix)
variable pix
Nvar borneHoughIntercepts1, intervalleHoughIntercepts, storageHoughResolution
return borneHoughIntercepts1 + pix * (intervalleHoughIntercepts/storageHoughResolution)
end

function pix2slo(pix)
variable pix
Nvar borneHoughSlopes1, intervalleHoughSlopes, storageHoughResolution
return borneHoughSlopes1 + pix * (intervalleHoughSlopes/storageHoughResolution)
end

function slo2pix(slo)
variable slo
Nvar borneHoughSlopes1, intervalleHoughSlopes, storageHoughResolution
return (slo-borneHoughSlopes1)*(storageHoughResolution/intervalleHoughSlopes)
end

function inter2pix(inter)
variable inter
Nvar borneHoughIntercepts1, intervalleHoughIntercepts, storageHoughResolution
return (inter-borneHoughIntercepts1)*(storageHoughResolution/intervalleHoughIntercepts)
end

Function changeHTreso(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		miseajourtransfo()
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function calcHTinterGaps(inter,slo)
variable inter,slo
wave listeMasseFacteursRosePentes
wave/T listeFacteursRosePentes
nvar resoenmasse,borneHoughSlopes2, borneHoughSlopes1, borneHoughIntercepts2, borneHoughIntercepts1, storageHoughResolution,intervalleHoughSlopes,intervalleHoughIntercepts,inco
nvar precslo,precinter
Duplicate/O listeMasseFacteursRosePentes HTinterGaps,HTsloGaps
HTinterGaps=inter2pix(inter-slo*listeMasseFacteursRosePentes+defmass(listeMasseFacteursRosePentes))
HTsloGaps=slo2pix(slo)
duplicate/O HTinterGaps indexgaps
Makeindex/R HTinterGaps indexgaps
variable n=numpnts(HTsloGaps),k
Make/O/D/N=(3*n) histohoughaffichagehoriz1, histohoughaffichagehoriz2
histohoughaffichagehoriz2=storageHoughResolution*(mod(x,3)!=0)
histohoughaffichagehoriz2/=(mod(x+1,3)!=0)
string letag=""
SetAxis/W=etudeTransformationHough#imageResultatHough left 0,storageHoughResolution
SetAxis/W=etudeTransformationHough#imageResultatHough bottom 0,storageHoughResolution
selecHT(inter,slo,precinter,precslo,inco)
If(winType("zoomht")!=0)
	SetAxis/W=zoomHT/R left wavemin(roiHTinter),wavemax(roiHTinter)
	SetAxis/W=zoomHT/R top wavemin(roiHTslo),wavemax(roiHTslo)
endIf
	if (stringmatch(TraceNameList("etudeTransformationHough#imageResultatHough", ";", 1 ),"*HTinterGaps*")==0)
		AppendToGraph/W=etudeTransformationHough#imageResultatHough HTinterGaps vs HTsloGaps
		ModifyGraph/W=etudeTransformationHough#imageResultatHough mode(HTinterGaps)=3,marker(HTinterGaps)=9
		ModifyGraph/W=etudeTransformationHough#imageResultatHough rgb(HTinterGaps)=(65280,0,0)
		AppendToGraph/W=etudeTransformationHough#imageResultatHough roiHTinter vs roiHTslo
	endif
for(k=0;k<n;k+=1)
	letag="letag"+num2str(k)
	Tag/W=etudeTransformationHough#imageResultatHough/C/N=$letag/L=3/X=(30*(1-2*(mod(k,2)==0)))/Y=(100*(round((n-k+2)/2)*storageHoughResolution/(n)-inter2pix(inter))/storageHoughResolution)/TL={lineRGB=(65280,0,0)}  HTinterGaps, indexgaps[k],listeFacteursRosePentes[indexgaps[k]]
endfor
histohoughaffichagehoriz1=inter2pix(	inco*(slo-pix2slo(histohoughaffichagehoriz2))+pix2inter(HTinterGaps[floor(x/3)])	)
end

function selecHT(inter,slo,preci,precs,m)
variable inter,slo,preci,precs,m
Make/O/D/N=5 roiSlo, roiInter, roiHTslo, roiHTinter
roiSlo={slo-precs,slo-precs,slo+precs,slo+precs,slo-precs}
roiInter={m*precs+inter+preci,m*precs+inter-preci,-m*precs+inter-preci,-m*precs+inter+preci,m*precs+inter+preci}
roiHTslo=slo2pix(roiSlo)
roiHTinter=inter2pix(roiInter)
end

function getHTroiMatch(inter,slo,preci,precs,m)
variable inter,slo,preci,precs,m
variable k,j,n,t=ticks, test=0
wave slopes, intercepts, modules,centremasse, wave0stripped, wave1stripped, dataindex
n=dimsize(slopes,0)
Make/O/D/N=0 lesk,lesj, roi0,roi1
for(k=0;k<n;k+=1)
	for(j=1+k;j<n;j+=1)
		if(slopes[j][k]>=(slo-precs))
			if(slopes[j][k]<=(slo+precs))
				if(mod(intercepts[j][k],1)>=mod((m*(slo-slopes[j][k])+inter-preci),1))
					if(mod(intercepts[j][k],1)<=mod((m*(slo-slopes[j][k])+inter+preci),1))
						test=1
						insertpoints 0,1, lesk, lesj
						lesk[0]=k
						lesj[0]=j
						//doupdate
						//j=n
					endif
				endif
			endif
		endif
	endfor
	if(test)
		insertpoints 0,1, roi0,roi1
		roi0[0]=wave0stripped[dataindex[k]]
		roi1[0]=wave1stripped[dataindex[k]]
		test=0
	endif
endfor
duplicate/O lesk, modtrans, pentrans, interceptrans, centremassetrans
modtrans=modules[lesj][lesk]
pentrans=slopes[lesj][lesk]
interceptrans=intercepts[lesj][lesk]
centremassetrans=centremasse[lesj][lesk]
sort modtrans, modtrans, pentrans, interceptrans, lesk,lesj, centremassetrans
duplicate/O roi1,roi2
sort roi0, roi0,roi1
roi2=defmass(roi0)
print ticks-t
end

function TraceRoiXpress(inter,slo,preci,precs,m)
variable inter,slo,preci,precs,m
variable k,j,n,t=ticks
wave slopes, intercepts, modules,centremasse, wave0stripped, wave1stripped, dataindex//,roi0,roi1
n=dimsize(slopes,0)
Make/O/D/N=0 roi0,roi1, roipnts
for(k=0;k<n;k+=1)
	for(j=1+k;j<n;j+=1)
		if(slopes[j][k]>=(slo-precs))
			if(slopes[j][k]<=(slo+precs))
				if(mod(intercepts[j][k],1)>=mod((m*(slo-slopes[j][k])+inter-preci),1))
					if(mod(intercepts[j][k],1)<=mod((m*(slo-slopes[j][k])+inter+preci),1))
						insertpoints 0,1, roi0,roi1,roipnts
						roipnts[0]=dataindex[k]
						roi0[0]=wave0stripped[dataindex[k]]
						roi1[0]=wave1stripped[dataindex[k]]
						j=n
					endif
				endif
			endif
		endif
	endfor
endfor
duplicate/O roi1,roi2
sort roi0, roi0,roi1
roi2=defmass(roi0)
print ticks-t
end

Function SetVarProc_6(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
variable/g treedepth
nvar inco
Variable dval = sva.dval
String sval = sva.sval
	treedepth=dval
	controlinfo/W=dmvm popup0
	generatetree(S_value,treedepth,inco)
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
		case -1: // control being killed
			break
	endswitch

	return 0
End

function HoughField(w,x,y) : fitfunc
wave w
variable x
variable y
return w[0]*(exp(-(y-(-w[1]*(x-w[2])+w[3]))^2/w[4]^2))/(1+(x-w[2])^2/w[5]^2)
end

function MultiHoughField(w,x,y) : fitfunc
wave w
variable x
variable y
variable res=0,k,n=11, gap=26.4
for(k=0;k<n;k+=1)
	res+= w[0]*(exp(-(y-(-w[1]*(x-w[2])+w[3]+(k-(n-1)/2)*gap))^2/w[4]^2))/(1+(x-w[2])^2/w[5]^2)
endfor
return res
end

Function ListPROPreturntoelaborator(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
wave/Z molecule
string lawave
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listwave = lba.listWave
	WAVE/Z selwave = lba.selWave
nvar inco,treedepth
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 3: // double click
			break
		case 4: // cell selection
			controlinfo/W=attributeur list0
			if(v_value>-1)
				lawave=stringfromlist(v_value,wavelist("molprop_*",";",""))
				molecule=0
				tuelesspec()
				duplicate/O $lawave molecule
				genestringmol()
				controlinfo/W=dmvm popup0
				generatetree(S_value,treedepth,inco)
			endif
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
	endswitch

	return 0
End

Function PopMenuProc_3(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
nvar treedepth, inco
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			generatetree(popStr,treedepth,inco)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function MolbagAddTo(molecule,charge,nom,cible,bag)
wave molecule,cible
variable charge
string nom, bag
variable t=ticks
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
variable n
if(waveexists($smolbag) && numpnts(molecule)==114 && numpnts(cible)==5)
	wave molbag=$smolbag
	wave chargebag=$schargebag
	wave/T listbag=$slistbag
	wave ciblebag=$sciblebag
	n=dimsize(molbag,1)
	if(n==0)
		make/O/D/N=(114,1) $smolbag
		make/O/D/N=(1,5) $sciblebag
	else
	insertpoints/M=1 n,1,molbag
	insertpoints/M=0 n,1,ciblebag
	endif
	molbag[][n]=molecule[x]
	insertpoints/M=0 n,1,chargebag,listbag
	chargebag[n]=charge
	listbag[n]=nom
	ciblebag[n][]=cible[y]
elseif(waveexists($smolbag)==0 && numpnts(molecule)==114 && numpnts(cible)==5)
	Make/O/D/N=(114,1) $smolbag
	Make/O/D/N=(1,5) $sciblebag
	Make/O/D/N=1 $schargebag
	Make/O/T/N=1 $slistbag
	wave molbag=$smolbag
	wave chargebag=$schargebag
	wave/T listbag=$slistbag
	wave ciblebag=$sciblebag
	molbag[][0]=molecule[x]
	ciblebag[0][]=cible[y]
	chargebag[0]=charge
	listbag[0]=nom
endif
//print "Temps mis pour ajouter une molécule au Molbag :", ticks-t
end

function MolbagDecant(from,to,which) //transvase
string from,to
wave which
variable n=(numpnts(which))
MolbagExtractor(from,which,to)
MolbagWithdrawFrom(which,from)
end

function MolbagWithdrawFrom(mol2kill,bag)
wave mol2kill// index of molecules to kill
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
if(waveexists(molbag))
variable n=numpnts(mol2kill),k
sort/R mol2kill, mol2kill
for(k=0;k<n;k+=1)
	deletepoints/M=1 mol2kill[k],1,molbag
	deletepoints/M=0 mol2kill[k],1,chargebag,listbag,ciblebag
endfor
endif
end

function MolbagEmpty(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
if(waveexists(molbag))
variable n=numpnts(chargebag)
deletepoints/M=1 0,n,molbag
deletepoints/M=0 0,n,chargebag,listbag,ciblebag
endif
end

function Molbag2generator(bag, item)
string bag
variable item
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
wave molecule
nvar charge
molecule=molbag[x][item]
charge=chargebag[item]
tuelesspec()
genestringmol()
end

function MolbagDetarget(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable k,n=dimsize(molbag,1)
for(k=0;k<n;k+=1)
	extract/O/FREE molbag,molecule,y==k
	ciblebag[k][1]=mol2massmono(molecule,chargebag[k])
	ciblebag[k][3]=0
endfor
end

Function PopMenuAddFromMolbags(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			addthemolbag(popstr)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_10(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave index_formules, list_targets
wave/T list_formules
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			edit/K=1 list_formules, index_formules, list_targets
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

/////////AdvancedMolManager
Function InitBagList(quant)
variable quant
string lasel="seladvmolman"+num2str(quant)
string lalist="advmolmanlist"+num2str(quant)
variable n=itemsinlist(wavelist("Molbag_*",";",""))
Make/O/T/N=6 advmolmantitles={"","Names","Points","Spec","DM",""}
Make/O/I/N=(n,6) $lasel=0
Make/O/T/N=(n,6) $lalist=""
wave seladvmolman=$lasel
wave/T advmolmanlist=$lalist
seladvmolman[][0]=64
seladvmolman[][1,2]=0
seladvmolman[][3,5]=32
advmolmanlist[][1]=stringfromlist(x,replacestring("Molbag_",wavelist("Molbag_*",";",""),""))
advmolmanlist[][2]=num2str(numpnts($stringfromlist(x,wavelist("Chargebag_*",";",""))))
end

function MajBagList(quant)
variable quant
string lasel="seladvmolman"+num2str(quant)
string lalist="advmolmanlist"+num2str(quant)
wave seladvmolman=$lasel
wave/T advmolmanlist=$lalist
//vérifier que ce qui est dans la liste est tout ce qui existe
//Priorité à ce qui existe, i.e. les membres de la liste uniquement sont delete
string lesbags=wavelist("Molbag_*",";","")
variable j=itemsinlist(lesbags), k, i, n=dimsize(seladvmolman,0),bot
make/O/D/N=(j,n) molbagmatch
molbagmatch=stringmatch(advmolmanlist[y][1],stringfromlist(x,replacestring("Molbag_",lesbags,"")))
MatrixOP/O testi=sumrows(molbagmatch)
MatrixOP/O testk=(sumcols(molbagmatch))^t
for(k=n-1;k>-1;k-=1)
	if(testk[k]==0)
		deletepoints/M=0 k,1,advmolmanlist,seladvmolman
	endif
endfor
for(i=0;i<j;i+=1)
	if(testi[i]==0)
		bot=dimsize(seladvmolman,0)
		insertpoints/M=0 bot,1, advmolmanlist,seladvmolman
		seladvmolman[bot][0]=64
		seladvmolman[bot][1,2]=0
		seladvmolman[bot][3,5]=32
		advmolmanlist[bot][1]=stringfromlist(i,replacestring("Molbag_",lesbags,""))
	endif
endfor
advmolmanlist[][2]=num2str(numpnts($("Chargebag_"+advmolmanlist[x][1])))
killwaves/Z molbagmatch//,testi,testk
end

function KillBagList(quant)
variable quant
string lasel="seladvmolman"+num2str(quant)
string lalist="advmolmanlist"+num2str(quant)
killwaves $lasel,$lalist
end

Function InitMolList(quant,bag)
variable quant
string bag
string lasel="selmol"+num2str(quant)
string lalist="mollist"+num2str(quant)
string latitle="moltitles"+num2str(quant)
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1)
Make/O/I/N=(n,11) $lasel=0
Make/O/T/N=(n,11) $lalist=""
Make/O/T/N=11 $latitle=""
wave selmol=$lasel
wave/T mollist=$lalist
wave/T moltitles=$latitle
moltitles={"","\f01"+bag,"Calc","Target","Bias (ppm)","Intensity","Cost Func.","Y","","",""}
selmol[][0]=64
selmol[][1]=0
selmol[][7]=0
selmol[][8,10]=32
mollist[][1]=listbag[x]
mollist[][3]=num2str(ciblebag[x][1])
mollist[][4]=num2str(ciblebag[x][3])
mollist[][2]=num2str(ciblebag[x][1]/(1-ciblebag[x][3]/1e6))
mollist[][5]=num2str(ciblebag[x][2])
mollist[][6]=num2str(ciblebag[x][0])
mollist[][7]=num2str(ciblebag[x][4])
end

function KillMolList(quant)
variable quant
string lasel="selmol"+num2str(quant)
string lalist="mollist"+num2str(quant)
string latitle="moltitles"+num2str(quant)
killwaves $lasel,$lalist,$latitle
end

Function PilotBagList(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
 variable quant=str2num(replacestring("advmolmanpanel",lba.win,""))
//print lba
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			majbaglist(quant)
			if(lba.eventMod==16 || lba.eventMod==17)
				if(row>-1 && col<3)
					string MenuStr=InitAdvMolManMenu(selWave,listWave)
					wave res
					PopupContextualMenu MenuStr
					if(V_flag==1)
						MolbagDBE2fifth(lba.listwave[row][1])
					// **************
					// ADD CW 12032019 for adding RT informations
					elseif(V_flag==2)
						// Open a window for picking the corresponding treated chromato
						// Check with the masses for matches
						// Add the RT to ciblebag[][4]
						string chromalist,sRTval
						Prompt chromalist,"Choose chromatogram for RT piking",popup,wavelist("*_11chroma",";","")
						DoPrompt/HELP="You have to select the chromatogram name where the system will extract the couple [mass;RT] and try to match it with the selected attribution" "RT selection", chromalist
						if (V_Flag)
							return -1	// User canceled
						else
							string sciblebag="ciblebag_"+lba.listwave[row][1]
							wave ciblebag=$sciblebag
							// Put Nan in Y_values
							ciblebag[][4]=Nan
							Duplicate/O/FREE ciblebag, transciblebag 
							extract/O/FREE transciblebag, transmass,y==1
							extract/O/FREE transciblebag, transerror,y==3
							//traitement de *_11chroma
							Duplicate/O/FREE $chromalist, target 
							extract/O/FREE target,transtargetmass,y==1
							make/O/U/FREE/N=(numpnts(transtargetmass)) targetmass
							targetmass=transtargetmass
							extract/O/FREE target,targetRT,y==3
							extract/O/FREE target,targetint,y==2
							variable j,l,pntattrib=numpnts(transmass),ciblemass, minciblemass, maxciblemass, error
							for(j=0;j<pntattrib;j+=1) //Bouclage sur les points de l'attrib
								sRTval=""
								ciblemass=transmass[j]
								error=transerror[j]
								minciblemass=transmass[j]-transmass[j]*abs(transerror[j])*10^(-6)
								maxciblemass=transmass[j]+transmass[j]*abs(transerror[j])*10^(-6)
								// Index des attrib ayant la masse demandée
								extract/O/FREE/INDX transmass,transcibleindx,transmass==ciblemass
								// Index des fit ayant la masse demandée
								extract/O/FREE/INDX targetmass,transtargetindx,targetmass<=maxciblemass && targetmass>=minciblemass
								// On vérifie le nombre de points trouvés dans le fit
								if(numpnts(transtargetindx)==0)
									// Do nothing
								elseif(numpnts(transtargetindx)==1)
									// Put in ciblebag[][4] the RT value with 2 significant digits + rewrite intensity value
									sprintf sRTval,"%.2f",targetRT[transtargetindx[0]]
									ciblebag[j][4]=str2num(sRTval)
									ciblebag[j][2]=targetint[transtargetindx[0]]
								elseif(numpnts(transtargetindx)>1)
									l=0
									for(l=0;l<numpnts(transtargetindx);l+=1)
										// Put in ciblebag[][4] the RT value with 2 significant digits + rewrite intensity value
										sprintf sRTval,"%.2f",targetRT[transtargetindx[l]]
										ciblebag[j+l][4]=str2num(sRTval)
										ciblebag[j+l][2]=targetint[transtargetindx[l]]
									endfor
									j=j+l-1
								else
									print "Echec - Manque une condition"
								endif
							endfor
						endif
					elseif(V_flag>2)
						//print V_flag, res
						MolbagElement2fifth(lba.listwave[row][1],res[V_flag-3])
					endif
					// END CW 12032019
					// **************
					   
					  
				elseif(row==-1 && (col==3 || col==4))
					MenuStr="Labelize Name;Name and Mass;Name, Mass and Error;No labels;"
					PopupContextualMenu MenuStr
					switch(V_flag)
						case 1:
							EtikTrace(col2win(col),"Traces_","",0,0,1)
							break
						case 2:
							EtikTrace(col2win(col),"Traces_","",1,0,1)
							break
						case 3:
							EtikTrace(col2win(col),"Traces_","",1,1,1)
							break
						case 4:
							EtikTrace(col2win(col),"Traces_","",0,0,0)
							break
					endswitch
				endif
			endif
			if(row>-1)
				InitMolList(quant,lba.listwave[row][1])
			endif
			break
		case 2:
			if(row==-1 && col<3)
				resortDATAbycol(selwave,listwave,col)
			elseif(row>-1 && col>2)
				OperateTraceOnGraph(listwave[row][1],col2win(col),(selwave[row][col]==48))
			endif
			break
		case 3: // double click
			if(row==-1 && col>2)
				selwave[][col]=32+16*(selwave[x][col]==32)
				MajTraceOnGraph(selwave,listwave)
			elseif(col==1)
				InitAdvMolManGraph(quant)
				MajAdvMolManGraph(quant)
			endif
			break
		case 4: // cell selection
			InitMolList(quant,lba.listwave[row][1])
			MajAdvMolManGraph(quant)
			if(lba.eventMod==9)//with ctrl
				if(col==1 || col==2)
					InitAdvMolManPanel(itemsinlist(winlist("AdvMolManPanel*",";","")),lba.listwave[row][1])
					autopositionwindow/M=1/R=$lba.win
				endif
			elseif(lba.eventMod==5)//with maj
				InitAdvMolManOper(quant)
			endif
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 12://keystroke
			if(row==127)
			extract listWave,killme,(y==1 && mod(selWave[x][0],2))
			variable k,n=dimsize(killme,0)
				for(k=0;k<n;k+=1)
					execute "MolbagKill(killme["+num2str(k)+"])"
				endfor
				majbaglist(quant)
				//InitMolList(quant,lba.listwave[v_value][1])
				//listbox baglist selrow=v_value-1
			endif
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			if(col==0)
				string item1="Fully decant "+listwave[row][1]+" to...;"
				string item2="Duplicate "+listwave[row][1]+";"
				string item3="Put "+listwave[row][1]+" to ROI;"
				string item4="Export "+listwave[row][1]+" table;"
				string item5="Update "+listwave[row][1]+" Mass biases and Intensities;"
				string item6="Operate to all items in "+listwave[row][1]+";"
				string item7="Rename;"
				string item8="Detarget"
				popupcontextualmenu item1+item2+item3+item4+item5+item6+item7+item8
				switch(V_flag)
					case 1:
						string newname
						prompt newname,"Decant "+listwave[row][1]+" into the following list",popup, replacestring("Molbag_",wavelist("Molbag_*",";",""),"")
						Doprompt "Chose a destination list",newname
						if(1-V_flag)
							MolbagDecant(listwave[row][1],newname,Molbag2FullIndx(listwave[row][1]))
							majbaglist(quant)
						endif
						break
					case 2:
						newname=listwave[row][1]+num2str(itemsinlist(wavelist("Molbag_"+listwave[row][1]+"*",";","")))
						prompt newname,"Give name below or overwrite existing:"
						Doprompt "Enter a name for a list",newname
						if(1-V_flag)
							MolbagDuplicate(listwave[row][1],newname)
							majbaglist(quant)
						endif
						break
					case 3:
						Molbag2ROI(listwave[row][1])
						break
					case 4:
						string stable=listwave[row][1]+"_table"
						Duplicate/O Molbag2Table(listwave[row][1]) $stable
						edit/K=1 titles
						edit/K=1 $stable
						break
					case 5:
						MolbagRescale(listwave[row][1])
						break
					case 6:
						string str
						prompt str, "The following line will be read, the stoichiometry added and the charge applied to the whole list"
						Doprompt "Enter a formula",str
						if(1-V_flag)
							MolbagOperateToALL(listwave[row][1],str2mol(str),str2charge(str))
						endif
						break
					case 7:
						newname=listwave[row][1]+num2str(itemsinlist(wavelist("Molbag_"+listwave[row][1]+"*",";","")))
						prompt newname,"Give name below or overwrite existing:"
						Doprompt "Enter a name for a list",newname
						if(1-V_flag)
							MolbagRename(listwave[row][1],newname)
						endif
						break
					case 8:
						MolbagDetarget(listwave[row][1])
						break
				endswitch
				selwave[row][col]=64
			endif
			MajAdvMolManGraph(quant)
			break
	endswitch
	return 0
End

Function PilotMolList(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
variable quant=str2num(replacestring("advmolmanpanel",lba.win,""))
string latitle="moltitles"+num2str(quant), name=""
string opername="AdvMolManOper"+num2str(quant)
wave/T moltitles=$latitle
//print lba
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
			if(row==-1 && col<8)
				name=moltitles[1]
				name=name[4,strlen(name)]
				//print name,col
				MolbagResort(name,col)
				InitMolList(quant,name)
				MajAdvMolManGraph(quant)
			elseif(row==-1 && col>7)
				selwave[][col]=32+16*(selwave[x][col]==32)
			endif
			break
		case 3: // double click
			if(col==1)
			wave molecule
			nvar charge
			string str=mol2instr(molecule,charge)
			prompt str,"Respect the format: ElementA#ElementB# +CHARGE, e.g. CH4 +1; C23H107N -12"
			doprompt "Edit the list item",str
			if(V_flag)
			else
				name=moltitles[1]
				name=name[4,strlen(name)]
				string smolbag="Molbag_"+name
				string schargebag="Chargebag_"+name
				string slistbag="Listbag_"+name
				string sciblebag="Ciblebag_"+name
				wave molbag=$smolbag
				wave chargebag=$schargebag
				wave/T listbag=$slistbag
				wave ciblebag=$sciblebag
				molbag[][row]=str2mol(str)[x]
				chargebag[row]=str2charge(str)
				listbag[row]=mol2str(str2mol(str),str2charge(str))
				Make/O/D/N=5/FREE transcible
				transcible={0,mol2massmono(str2mol(str),str2charge(str)),ciblebag[row][2],0,0}
				ciblebag[row][]=transcible[y]
				InitMolList(quant,name)
			endif
			endif
			break
		case 4: // cell selection
			name=moltitles[1]
			name=name[4,strlen(name)]
			Molbag2generator(name, row)
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 12://keystroke
			if(row==127)
				extract/FREE/INDX selwave,killme,selWave[x][y]>64 && y==0
				name=moltitles[1]
				name=name[4,strlen(name)]
				MolbagWithdrawFrom(killme,name)
				InitMolList(quant,name)
				MajBagList(quant)
			endif
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch

	return 0
End

Function AdvMolManGraphPopMenuProc_4(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
 variable quant=str2num(replacestring("advmolmangraph",pa.win,""))
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			MajAdvMolManGraph(quant)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


function InitAdvMolManPanel(quant,bag)
variable quant
string bag
InitBagList(quant)
InitMolList(quant,bag)
string panelname="AdvMolManPanel"+num2str(quant)
string molsel="selmol"+num2str(quant)
string mollist="mollist"+num2str(quant)
string bagsel="seladvmolman"+num2str(quant)
string baglist="advmolmanlist"+num2str(quant)
string latitle="moltitles"+num2str(quant)
wave selWave=$bagsel
Newpanel/K=1/N=$panelname/W=(10,10,410,410)
ListBox baglist,pos={1,1},size={400,200},proc=PilotBagList
ListBox baglist,listWave=$baglist,selWave=$bagsel
FindValue/TEXT=bag/TXOP=4 advmolmanlist0
ListBox baglist,titleWave=root:advmolmantitles,mode= 4//,selRow= v_value-floor(v_value/dimsize($baglist,0))*dimsize($baglist,0)
selWave[v_value-floor(v_value/dimsize($baglist,0))*dimsize($baglist,0)][0]+=1
ListBox baglist,widths={15,150,50,15,15,15}
ListBox mollist,pos={1,200},size={400,200},labelBack=(47872,47872,47872),proc=PilotMolList//,fSize=10
ListBox mollist,frame=0,listWave=$mollist,selWave=$molsel
ListBox mollist,titleWave=$latitle,mode= 9,special= {0,20,0}
ListBox mollist,widths={0,75,0,30,40,40,0,20,0,0,15},userColumnResize= 1
setwindow $panelname, hook(hAdvMolMan)=AdvMolManPanelHook
end

function MajAdvMolManPanel(quant)
variable quant
string panelname="AdvMolManPanel"+num2str(quant)
getwindow $panelname wsizeDC
variable width=v_right,height=v_bottom
string mollist="mollist"+num2str(quant)
string baglist="advmolmanlist"+num2str(quant)
ListBox baglist,pos={1,1}, win=$panelname, size={width,min(200,height/2)}
ListBox mollist,pos={1,min(200,height/2)}, win=$panelname, size={width,height-min(200,height/2)}
end

Function AdvMolManPanelHook(s)
	STRUCT WMWinHookStruct &s
string graphname=replacestring("advmolmanpanel",s.winName,"AdvMolManGraph")
string panelname=s.winName
variable quant=str2num(replacestring("advmolmanpanel",s.winName,""))
	Variable hookResult = 0

	switch(s.eventCode)
		case 2:
			killmollist(quant)
			killbaglist(quant)
			break
		case 6:
			MajAdvMolManPanel(quant)
			break
		case 12:
			if(wintype(graphname))
				autopositionwindow/M=0/R=$panelname $graphname
			endif
			break
	endswitch
	return hookResult		// 0 if nothing done, else 1
End

Function MajTagAdvMolMan(quant,name,hitpoint)
variable quant, hitpoint
string name
string opername="AdvMolManOper"+num2str(quant)
String graphname="AdvMolManGraph"+num2str(quant)
variable k,l,m,n
string lafenetre, lestraces,latrace, lacouleur
for(k=0;k<2;k+=1)
	lafenetre=stringfromlist(k,opername+";"+graphname+";")
	lestraces=tracenamelist(lafenetre,";",1)
	n=itemsinlist(lestraces)
	for(l=n-1;l>-1;l-=1)
		latrace=stringfromlist(l,lestraces)
		if(stringmatch(trace2bagname(latrace),name)==1)
			string slistbag="Listbag_"+name
			wave/T listbag=$slistbag
			Tag/W=$lafenetre/C/N=text0/A=LB/Y=(5)/X=(5)/L=2 $latrace , hitpoint, listbag[hitpoint]
		endif
	endfor
endfor
end

Function AdvMolManGraphHook(s)
	STRUCT WMWinHookStruct &s
string panelname=replacestring("AdvMolManGraph",s.winName,"advmolmanpanel")
string graphname=s.winName
variable quant=str2num(replacestring("AdvMolManGraph",s.winName,""))
string molsel="selmol"+num2str(quant)
wave selmol=$molsel
variable hitpoint=str2num(stringbykey("HITPOINT", TraceFromPixel(s.mouseLoc.h, s.mouseLoc.v,""))),k
string tracename=stringbykey("TRACE", TraceFromPixel(s.mouseLoc.h, s.mouseLoc.v,""))
string infotrace=traceinfo(graphname,tracename,0)
string name=trace2bagname(tracename)
string tracelist=tracenamelist(graphname,";",1)
string hittracelist=""
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
string sClass="Clu_Cla9_"+name+"Oper"+num2str(quant)
if(waveexists($sClass))
	wave Class=$sClass
endif
if(strlen(name))
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
endif
string popupMenuC=""
variable nbclass=3
for(k=0;k<itemsinlist(tracelist);k+=1)
	if(strsearch(stringfromlist(k,tracelist),name,0)!=-1)
		hittracelist+=stringfromlist(k,tracelist)+";"
	endif
endfor
Variable hookResult = 0
Make/FREE/T/N=10 cursname={"A","B","C","D","E","F","G","H","I","J"}
for(k=k;k<10;k+=1)
	Cursor/K/W=$graphname $cursname[k]
endfor
if(strlen(name))
	//Tag/W=$graphname/C/N=text0/A=LB/Y=(5)/X=(5)/L=2 $tracename , hitpoint, listbag[hitpoint]+" "+num2str(ciblebag[hitpoint][3])+name
	MajTagAdvMolMan(quant,name,hitpoint)
endif
for(k=0;k<itemsinlist(hittracelist);k+=1)
	Cursor/W=$graphname/A=0/H=0/S=0/C=(26368,0,52224) $cursname[k], $stringfromlist(k,hittracelist), hitpoint
endfor
Cursor/W=$graphname/M/A=0/H=2/S=0/L=1 A
	switch(s.eventCode)
		case 2:
			string Clustering2kill=wavelist("Clu_*Oper"+num2str(quant),";",""), wave2kill
			variable n2kill
			n2kill=itemsinlist(Clustering2kill)
			for(k=0;k<n2kill;k+=1)
				wave2kill=stringfromlist(k,Clustering2kill)
				killwaves/Z $wave2kill
			endfor
			break
		case 3:
		selmol[][0]=64
		selmol[hitpoint][0]=65
		listbox mollist win=$panelname, row=max(-1,hitpoint-2)
		if(s.eventmod==16 ||s.eventmod==17 && stringmatch(tracename,"")==0 && strsearch(tracename,"Ciblebag_",0)!=-1 && stringmatch(stringbykey("YRANGE",infotrace),"[*][3]"))
			popupmenuc="Cluster Mass Bias distance to polynomials for "+name+";"+"Calculate "+num2str(Ciblebag[hitpoint][1])+";"
			if(1)
				popupmenuc+="Delete this class;"
			endif
			PopupContextualMenu popupmenuc
			switch(V_flag)	// numeric switch
				case 1:
					ClassBiasPoly(name,6,nbclass,quant)//name,Order,nbClust,quant
					wave Class=$sClass
					ModifyGraph/W=$graphname zColor($tracename)={Class,-1,*,Classification,0}, marker($tracename)=19, useMrkStrokeRGB($tracename)=1//, zmrksize($tracename)={vorosize,-1,*,wavemin(vorosize),wavemax(vorosize)}
					break
				case 2:
					Calculate(Ciblebag[hitpoint][1],Ciblebag[hitpoint][2],Chargebag[hitpoint])
					Majbaglist(quant)
					break
				case 3:
					RemoveClassif(name,class,class[hitpoint])
					Majbaglist(quant)
					break
			endswitch
		endif
			break
		case 6:
			break
		case 12:
			break
		case 22:
		print "TBD"
			break
	endswitch
	return hookResult		// 0 if nothing done, else 1
End

function/S trace2bagname(tracename)
string tracename
tracename=replacestring("Ciblebag_",tracename,"")
tracename=replacestring("Molbag_",tracename,"")
tracename=replacestring("Ratio_",tracename,"")
tracename=replacestring("KnariesPos_",tracename,"")
if(strsearch(tracename,"#",0)!=-1)
	tracename=tracename[0,strsearch(tracename,"#",0)-1]
endif
if(strsearch(tracename,"Oper",0)!=-1)
	tracename=tracename[0,strsearch(tracename,"Oper",0)-1]
endif
return tracename
end

function InitAdvMolManGraph(quant)
variable quant
string graphname="AdvMolManGraph"+num2str(quant)
if(wintype(graphname)==0)
string panelname="AdvMolManPanel"+num2str(quant)
getwindow $panelname wsize
string molsel="selmol"+num2str(quant)
string mollist="mollist"+num2str(quant)
string bagsel="seladvmolman"+num2str(quant)
string baglist="advmolmanlist"+num2str(quant)
Display/K=1/N=$graphname/W=(0,0,607.5,v_bottom-v_top)
autopositionwindow/M=0/R=$panelname $graphname
ControlBar/W=$graphname 50
PopupMenu TopMenu mode=4,pos={1,1},value="Stoichiometry;Target Mass;Intensity;Relative mass bias;Y=f(X);Cost function",proc=AdvMolManGraphPopMenuProc_4
PopupMenu BottomMenu mode=1,pos={1,26},value="Stoichiometry;Target Mass;Intensity;Relative mass bias;Y=f(X);Cost function",proc=AdvMolManGraphPopMenuProc_4
PopupMenu AbsMenu title="vs", mode=1,pos={150,13},value="Target Mass;Intensity;Relative mass bias;Y=f(X);Cost function",proc=AdvMolManGraphPopMenuProc_4
Button OpenOper title="Open operation window", pos={400,13}, size={125,20}, proc=AdvMolManGraphButtonOpenOper
setwindow $graphname, hook(hAdvMolManG)=AdvMolManGraphHook
endif
end

Function AdvMolManGraphButtonOpenOper(ba) : ButtonControl
	STRUCT WMButtonAction &ba
variable quant=str2num(replacestring("AdvMolManGraph",ba.win,""))
	switch( ba.eventCode )
		case 2: // mouse up
			InitAdvMolManOper(quant)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function InitAdvMolManOper(quant)
variable quant
string panelname="AdvMolManPanel"+num2str(quant)
string opername="AdvMolManOper"+num2str(quant)
string graphname="AdvMolManGraph"+num2str(quant)
variable lasttab=0
if(wintype(panelname)==0 || wintype(opername)==1)
	controlinfo/W=$opername tab0
	lasttab=V_value
	dowindow/K $opername
	//return 0
endif
getwindow $panelname wsize
variable width=v_right-v_left, height=v_bottom-v_top, barsize=150
Display/K=1/N=$opername/W=(0,0,width,height)
autopositionwindow/M=0/R=$panelname $opername
if(wintype(graphname))
	autopositionwindow/M=0/R=$graphname $opername
endif
ControlBar/W=$opername barsize
setwindow $opername, hook(hAdvMolManO)=AdvMolManOperHook
getwindow $panelname wsizeDC
width=v_right-v_left
height=v_bottom-v_top
TabControl tab0, win=$opername, pos={1,1},size={width,barsize},proc=Pilottab, value=lasttab
TabControl tab0, win=$opername, tabLabel(0)="Clustering",tabLabel(1)="Ratio Plot",tabLabel(2)="K-naries",tabLabel(3)="Partitioning"
InitOperLists(quant)
setcontrolsOPER(lasttab,opername)
end

function setcontrolsOPER(n,win)
variable n
string win
variable quant=str2num(replacestring("AdvMolManOper",win,""))
//InitOperLists(quant)
string sListEltOper="ListEltOper"+num2str(quant)
string sSelEltOper="SelEltOper"+num2str(quant)
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sSelRatioOper="SelRatioOper"+num2str(quant)
string sListRatioOper="ListRatioOper"+num2str(quant)
string sClusteringTitlesOper="ClusteringTitlesOper"+num2str(quant)
string sSelClusteringOper="SelClusteringOper"+num2str(quant)
string sListClusteringOper="ListClusteringOper"+num2str(quant)
string sSelAlgoOper="SelAlgoOper"+num2str(quant)
string sListAlgoOper="ListAlgoOper"+num2str(quant)
string sSetVarNbclassOper="SetVarNbclassOper"+num2str(quant)
	switch(n)
		case 0:
		controlinfo/W=$win tab0
		SetVariable TabCon_Clustering_SetClus, win=$win, pos={v_left+170,v_top+20},size={v_width/4,v_height-25}, title="Classes",limits={1,inf,1},value=$sSetVarNbclassOper,frame=0,proc=ClusteringSetVar
		Listbox TabCon_Clustering_Metrics, win=$win, pos={v_left+3,v_top+20},size={160,v_height-25},widths={1,4,1,1},userColumnResize= 0, listWave=$sListClusteringOper, selWave= $sSelClusteringOper, titleWave=$sClusteringTitlesOper,frame=0,proc=ClusteringListBox
		Listbox TabCon_Clustering_Algo, win=$win, pos={v_left+170,v_top+40},size={v_width/4,v_height-45}, userColumnResize= 0, listWave=$sListAlgoOper, selWave= $sSelAlgoOper,frame=0,proc=ClusteringListBox
		//listbox TabCon_Clustering_Bagz, win=$win, pos={v_left+170+v_width/4,v_top+40},size={v_width/4,v_height-45}, listWave=$sListBagOper, selWave= $sSelBagOper,frame=0
		InitClusteringLists(quant)
		PopupMenu TabCon_Clustering_Color,  title="Color by",win=$win, pos={v_left+170+v_width/4,v_top+20}, mode=1,value=BuildColorListOper(), proc=ClusteringPopUp
			break
		case 1:
		controlinfo/W=$win tab0
		Listbox TabCon_Ratio_Formula, win=$win, pos={v_left+3,v_top+20},size={v_width-7,40},widths={50,150,25,25},userColumnResize= 0, listWave=$sListRatioOper, selWave= $sSelRatioOper,frame=0,proc=RatioListBox
		listbox TabCon_Ratio_Bagz, win=$win, pos={v_left+3,v_top+60},size={v_width-7,v_height-65}, listWave=$sListBagOper, selWave= $sSelBagOper,frame=0,proc=RatioListBox
		InitRatioLists(quant)
			break
		case 2:
		controlinfo/W=$win tab0
		listbox TabCon_Knaries_Elt, win=$win, pos={v_left+3,v_top+20},size={60,v_height-25}, listWave=$sListEltOper, selWave= $sSelEltOper, proc=KnariesListBox,frame=0
		listbox TabCon_Knaries_Bagz, win=$win, pos={v_left+63,v_top+20},size={v_width-68,v_height-25}, listWave=$sListBagOper, selWave= $sSelBagOper, proc=KnariesListBox,frame=0
		InitKnariesLists(quant)
			break
		case 3:
			break
	endswitch
end

Function ClusteringPopUp(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
 variable quant=str2num(replacestring("advmolmanoper",pa.win,""))
 string opername="AdvMolManOper"+num2str(quant)
 String graphname="AdvMolManGraph"+num2str(quant)
 variable k,l,m,n
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string lafenetre, lestraces,latrace, name, lacouleur
			//print claName2claNum(popStr)
			//Pour chaque fenetre et pour chaque trace, former le nom de la wave de classif
			for(k=0;k<2;k+=1)
				lafenetre=stringfromlist(k,opername+";"+graphname+";")
				lestraces=tracenamelist(lafenetre,";",1)
				n=itemsinlist(lestraces)
				for(l=0;l<n;l+=1)
					latrace=stringfromlist(l,lestraces)
					name=trace2bagname(latrace)
					lacouleur="Clu_Cla"+claName2claNum(popStr)+"_"+name+"Oper"+num2str(quant)
					Modifygraph/W=$lafenetre zColor($latrace)={$lacouleur,*,*,Classification}
				endfor
			endfor
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

function/S claName2claNum(claName)
string claName
string claNum=""
	strswitch(claName)
		case "Lloyd":
			claNum="0"
			break
		case "Kmeans":
			claNum="1"
			break
		case "FPnt":
			claNum="2"
			break
		case "Spectral":
			claNum="3"
			break
		case "Mass bias to polynome":
			claNum="9"
			break
	endswitch
return claNum
end 

Function ClusteringSetVar(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
 variable quant=str2num(replacestring("advmolmanoper",sva.win,""))
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			InitClusteringLists(quant)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ClusteringListBox(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
 variable quant=str2num(replacestring("advmolmanoper",lba.win,""))
//print lba
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
			break
		case 3: // double click
			break
		case 4: // cell selection
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			if(strlen(listWave[row][col])==0)
				deletepoints/M=0 row,1, listWave, selWave
			endif
			InitClusteringLists(quant)
			break
		case 12://keystroke
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			if(selwave[row][col]==80)
				Insertpoints/M=0 row,1, listWave, selWave
				selWave[row][]={{32},{2},{2},{2}}
				ListWave[row][]={{""},{"Mass"},{"1"},{"1"}}
				selWave[row+1][]={{64},{0},{0},{0}}
				ListWave[row+1][]={{""},{""},{""},{""}}
			else
				InitClusteringLists(quant)
			endif
			break
	endswitch
	return 0
End

Function KnariesListBox(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
 variable quant=str2num(replacestring("advmolmanoper",lba.win,""))
//print lba
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
			break
		case 3: // double click
			break
		case 4: // cell selection
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 12://keystroke
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
		InitKnariesLists(quant)
			break
	endswitch
	return 0
End

Function RatioListBox(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
 variable quant=str2num(replacestring("advmolmanoper",lba.win,""))
//print lba
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
			break
		case 3: // double click
			break
		case 4: // cell selection
			if(col==2)
				listwave[row][1]="1+C+N/2-H/2"
			elseif(col==3)
				listwave[row][1]="Mass"
			elseif(col==4)
				listwave[row][1]="C-(N+H)/2"
			elseif(col==5)
				listwave[row][1]="(H-3N)/2"
			endif
			InitRatioLists(quant)
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
		InitRatioLists(quant)
			break
		case 12://keystroke
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
		InitRatioLists(quant)
			break
	endswitch
	return 0
End

function/s wave2list(ws)
wave/T ws
string res=""
variable k,n=numpnts(ws)
for(k=0;k<n;k+=1)
	res+=ws[k]+";"
endfor
return res
end

function/S BuildColorListOper()
//extraire les numéros des classifications
string lesclass=wavelist("Clu_Cla*Oper*",";",""),laclass
variable k, n=itemsinlist(lesclass)
string res=""
string pseudores=""
string elt1
for(k=0;k<n;k+=1)
	laclass=stringfromlist(k,lesclass)
	splitstring/E="Clu_Cla([[:digit:]]+)" laclass, elt1
	if(strlen(listmatch(pseudores,elt1))==0)
		pseudores+=elt1+";"
		res+=claNum2claName(elt1)+";"
	endif
endfor
return res
end

function/S claNum2claName(claNum)
string claNum
string claName=""
	strswitch(claNum)
		case "0":
			claName="Lloyd"
			break
		case "1":
			claName="Kmeans"
			break
		case "2":
			claName="FPnt"
			break
		case "3":
			claName="Spectral"
			break
		case "9":
			claName="Mass bias to polynome"
			break
	endswitch
return claName
end

function InitOperLists(quant)
variable quant
string bagsel="seladvmolman"+num2str(quant)
string baglist="advmolmanlist"+num2str(quant)
string sIndxEltOper="IndxEltOper"+num2str(quant)
string sEltOper="EltOper"+num2str(quant)
string sSelEltOper="SelEltOper"+num2str(quant)
string sListEltOper="ListEltOper"+num2str(quant)
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sKnariesPoles="KnariesPolesOper"+num2str(quant)
string sSelRatioOper="SelRatioOper"+num2str(quant)
string sListRatioOper="ListRatioOper"+num2str(quant)
string sClusteringTitlesOper="ClusteringTitlesOper"+num2str(quant)
string sSelClusteringOper="SelClusteringOper"+num2str(quant)
string sListClusteringOper="ListClusteringOper"+num2str(quant)
string sSelAlgoOper="SelAlgoOper"+num2str(quant)
string sListAlgoOper="ListAlgoOper"+num2str(quant)
string sSetVarNbclassOper="SetVarNbclassOper"+num2str(quant)
string bagz=Ticked2List($bagsel,$baglist)
variable k,nbbag=itemsinlist(bagz)
Duplicate/O NonZeroElement(bagz), $sIndxEltOper, $sSelEltOper
Make/O/T/N=(dimsize($sIndxEltOper,0)) $sEltOper, $sListEltOper
make/O/D/N=(nbbag) $sSelBagOper
make/O/T/N=(nbbag) $sListBagOper
Make/O/D/N=(0,2) $sKnariesPoles
Make/O/D/N=(2,6) $sSelRatioOper
Make/O/T/N=(2,6) $sListRatioOper
Make/O/T/N=(4) $sClusteringTitlesOper
Make/O/D/N=(7+dimsize($sIndxEltOper,0),4) $sSelClusteringOper
Make/O/T/N=(7+dimsize($sIndxEltOper,0),4) $sListClusteringOper
Make/O/D/N=(3) $sSelAlgoOper
Make/O/T/N=(3) $sListAlgoOper
Make/O/D/N=1 $sSetVarNbclassOper=2
wave/T ListClusteringOper=$sListClusteringOper
wave SelClusteringOper=$sSelClusteringOper
wave/T ClusteringTitlesOper=$sClusteringTitlesOper
wave/T ListRatioOper=$sListRatioOper
wave SelRatioOper=$sSelRatioOper
wave IndxEltOper=$sIndxEltOper
wave/T EltOper=$sEltOper
wave/T elements
wave SelEltOper=$sSelEltOper
wave/T ListEltOper=$sListEltOper
wave SelBagOper=$sSelBagOper
wave/T ListBagOper=$sListBagOper
wave SelAlgoOper=$sSelAlgoOper
wave/T ListAlgoOper=$sListAlgoOper
ClusteringTitlesOper={"","Metrics","Weight","Norm"}
EltOper=elements[IndxEltOper[x]]
SelEltOper=48
ListEltOper=elements[IndxEltOper[x]]
ListBagOper=stringfromlist(x,bagz)
SelBagOper=48
SelRatioOper[][1]=2
ListRatioOper[][0]={"RatioY","RatioX"}
ListRatioOper[][1]={"1+C+N/2-H/2","Mass"}
ListRatioOper[][2]="DBE"
ListRatioOper[][3]="Mass"
ListRatioOper[][4]="Zubarev"
ListRatioOper[][5]="Meth.Excess"
ListClusteringOper[][1]=EltOper[x]
ListClusteringOper[dimsize($sIndxEltOper,0)][1]={"RatioY","RatioX","K-nariesY","K-nariesX","Error","Intensity"}
ListClusteringOper[][2,3]="1"
SelClusteringOper[][0]=32
SelClusteringOper[0,dimsize($sIndxEltOper,0)-1][0]+=16
SelClusteringOper[][2,3]=2
SelAlgoOper=48
ListAlgoOper={"Lloyd","Kmeans","FPnt"}
ListClusteringOper[dimsize(ListClusteringOper,0)-1][]=""
SelClusteringOper[dimsize(SelClusteringOper,0)-1][*]={{64},{2},{0},{0}}
end

function CleanupOper(quant)
variable quant
string opername="AdvMolManOper"+num2str(quant)
razTraces(opername)
string sIndxEltOper="IndxEltOper"+num2str(quant)
string sEltOper="EltOper"+num2str(quant)
string sSelEltOper="SelEltOper"+num2str(quant)
string sListEltOper="ListEltOper"+num2str(quant)
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sKnariesPoles="KnariesPolesOper"+num2str(quant)
string sSelRatioOper="SelRatioOper"+num2str(quant)
string sListRatioOper="ListRatioOper"+num2str(quant)
string sClusteringTitlesOper="ClusteringTitlesOper"+num2str(quant)
string sSelClusteringOper="SelClusteringOper"+num2str(quant)
string sListClusteringOper="ListClusteringOper"+num2str(quant)
string sSelAlgoOper="SelAlgoOper"+num2str(quant)
string sListAlgoOper="ListAlgoOper"+num2str(quant)
string sSetVarNbclassOper="SetVarNbclassOper"+num2str(quant)
Killwaves $sIndxEltOper, $sEltOper, $sSelEltOper, $sListEltOper, $sSelBagOper, $sListBagOper, $sKnariesPoles, $sSelRatioOper, $sListRatioOper,$sClusteringTitlesOper
Killwaves $sListClusteringOper, $sSelClusteringOper, $sSelAlgoOper, $sListAlgoOper, $sSetVarNbclassOper
string Knaries2kill=wavelist("KnariesPos_*Oper"+num2str(quant),";","")
variable k, n=itemsinlist(Knaries2kill)
string wave2kill
for(k=0;k<n;k+=1)
	wave2kill=stringfromlist(k,Knaries2kill)
	killwaves $wave2kill
endfor
string Ratio2kill=wavelist("Ratio_*Oper"+num2str(quant),";","")
n=itemsinlist(Ratio2kill)
for(k=0;k<n;k+=1)
	wave2kill=stringfromlist(k,Ratio2kill)
	killwaves $wave2kill
endfor
string Obs2kill=wavelist("Obs_*Oper"+num2str(quant),";","")
n=itemsinlist(Obs2kill)
for(k=0;k<n;k+=1)
	wave2kill=stringfromlist(k,Obs2kill)
	killwaves $wave2kill
endfor
string Clustering2kill=wavelist("Clu_*Oper"+num2str(quant),";","")
n=itemsinlist(Clustering2kill)
for(k=0;k<n;k+=1)
	wave2kill=stringfromlist(k,Clustering2kill)
	killwaves/Z $wave2kill
endfor
end

function InitClusteringLists(quant)
variable quant
string opername="AdvMolManOper"+num2str(quant)
if(wintype(opername)==0)
	return 0
endif
//razTraces(opername)
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sSelClusteringOper="SelClusteringOper"+num2str(quant)
string sListClusteringOper="ListClusteringOper"+num2str(quant)
string sSelAlgoOper="SelAlgoOper"+num2str(quant)
string sListAlgoOper="ListAlgoOper"+num2str(quant)
string sSetVarNbclassOper="SetVarNbclassOper"+num2str(quant)
wave SetVarNbclassOper=$sSetVarNbclassOper
wave/T ListBagOper=$sListBagOper
wave SelBagOper=$sSelBagOper
wave/T ListClusteringOper=$sListClusteringOper
wave SelClusteringOper=$sSelClusteringOper
wave/T ListAlgoOper=$sListAlgoOper
wave SelAlgoOper=$sSelAlgoOper
variable k,l,m,j
variable nbbag=dimsize(SelBagOper,0)
variable nbalgo=dimsize(SelAlgoOper,0)
variable nbmet=dimsize(SelClusteringOper,0)
variable nbmol=0
variable nbclass=0
variable metticked=0
Extract/FREE SelClusteringOper, tocount, SelClusteringOper==48
metticked=numpnts(tocount)
//Constituer les obs
string sObs="",sMolbag="",strans=""
wave/T elements
string selt=wave2list(elements)
variable scale
for(k=0;k<nbbag;k+=1)
	if(SelBagOper[k]==48)
	sObs="Obs_"+ListBagOper[k]+"Oper"+num2str(quant)
	sMolbag="Molbag_"+ListBagOper[k]
	nbmol=dimsize($sMolbag,1)
	Make/O/D/N=(nbmol,metticked) $sObs
	Make/FREE/O/D/N=(nbmol) FutureObs
	wave Obs=$sObs
	j=0
		for(m=0;m<nbmet;m+=1)
			if(SelClusteringOper[m][0]==48)
				strswitch(ListClusteringOper[m][1])
					case "RatioY":
					strans="Ratio_"+ListBagOper[k]+"Oper"+num2str(quant)
					wave trans=$strans
					FutureObs=trans[x][0]
						break
					case "RatioX":
					strans="Ratio_"+ListBagOper[k]+"Oper"+num2str(quant)
					wave trans=$strans
					FutureObs=trans[x][1]
						break
					case "K-nariesY":
					strans="KnariesPos_"+ListBagOper[k]+"Oper"+num2str(quant)
					wave trans=$strans
					FutureObs=trans[x][1]
						break
					case "K-nariesX":
					strans="KnariesPos_"+ListBagOper[k]+"Oper"+num2str(quant)
					wave trans=$strans
					FutureObs=trans[x][0]
						break
					case "Error":
					strans="Ciblebag_"+ListBagOper[k]
					wave trans=$strans
					FutureObs=trans[x][3]
						break
					case "Intensity":
					strans="Ciblebag_"+ListBagOper[k]
					wave trans=$strans
					FutureObs=trans[x][2]
						break
					default:
					variable WhichElt=WhichListItem(ListClusteringOper[m][1],selt,";" ,0,1)
					if(WhichElt==-1)
						Duplicate/O FutureObs ExpediantObs
						wave ExpediantObs
						string str2exec="ExpediantObs="+parse2eval(ListClusteringOper[m][1],ListBagOper[k])
						Execute str2exec
						FutureObs=ExpediantObs
						Killwaves ExpediantObs
					else
						strans=sMolbag
						wave trans=$strans
						FutureObs=trans[WhichElt][x]
					endif
						break
				endswitch
				scale=str2num(ListClusteringOper[m][2])/(wavemax(FutureObs)-wavemin(FutureObs))
				if(numtype(scale)==0)
					FutureObs*=scale
				endif
				Obs[][j]=FutureObs[x]
				j+=1
			endif
		endfor
		nbClass=SetVarNbclassOper[0]
		string Clustering2kill=wavelist("Clu_*Oper"+num2str(quant),";","")
		string wave2kill
		variable nclu2kill
		nclu2kill=itemsinlist(Clustering2kill)
		for(k=0;k<nclu2kill;k+=1)
			wave2kill=stringfromlist(k,Clustering2kill)
			killwaves $wave2kill
		endfor
		for(l=0;l<nbalgo;l+=1)
			if(SelAlgoOper[l]==48)
				string sClass="Clu_Cla"+num2str(l)+"_"+ListBagOper[k]+"Oper"+num2str(quant)
				string sCard="Clu_Car"+num2str(l)+"_"+ListBagOper[k]+"Oper"+num2str(quant)
				string sCenter="Clu_Cen"+num2str(l)+"_"+ListBagOper[k]+"Oper"+num2str(quant)
				strswitch(ListAlgoOper[l])
					case "Lloyd":
						Cluster(nbclass,Obs,2)
						Duplicate/O Means $sCenter
						Duplicate/O Card $sCard
						Duplicate/O Voronass $sClass
						break
					case "Kmeans":
						matrixOP/FREE toanal=Obs^t
						make/FREE/N=(metticked,nbClass) inw=enoise(1)
						Kmeans/CAN/INW=inw/OUT=2/TER=2/TERN=0 toanal
						wave M_KMClasses,W_KMMembers
						Make/O/D/N=(dimsize(M_KMClasses,1),dimsize(M_KMClasses,0)-1) $sCenter
						wave Center=$sCenter
						Center=M_KMClasses[y][x]
						Make/O/D/N=(dimsize(M_KMClasses,1)) $sCard
						wave Card=$sCard
						Card=M_KMClasses[dimsize(M_KMClasses,0)-1][x]
						Duplicate/O W_KMMembers $sClass
						break
					case "FPnt":
						FPClustering/NOR/MAXC=(nbClass) Obs
						wave W_FPClusterIndex, W_FPCenterIndex
						Duplicate/O W_FPClusterIndex $sClass
						wave Class=$sClass
						Make/O/D/N=(dimsize(W_FPCenterIndex,0),metticked) $sCenter
						wave Center=$sCenter
						Center=Obs[W_FPCenterIndex[x]][y]
						Make/O/D/N=(dimsize(W_FPCenterIndex,0)) $sCard
						wave Card=$sCard
						variable cfpnt
						for(cfpnt=0;cfpnt<dimsize(W_FPCenterIndex,0);cfpnt+=1)
							Extract/FREE W_FPClusterIndex,transitrans, W_FPClusterIndex==cfpnt
							Card[cfpnt]=sum(transitrans)
						endfor
						break
				endswitch
			endif
		endfor
	endif
endfor
end

function InitRatioLists(quant)
variable quant
string opername="AdvMolManOper"+num2str(quant)
if(wintype(opername)==0)
	return 0
endif
razTraces(opername)
string sRatio, sMolbag
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sListRatioOper="ListRatioOper"+num2str(quant)
wave/T ListBagOper=$sListBagOper
wave SelBagOper=$sSelBagOper
wave/T ListRatioOper=$sListRatioOper
variable nbbag=dimsize($sSelBagOper,0),nbmol,k,j,m,kpar,npar,V_flag=1
MAKE/FREE/O/D/N=10 red={65280,65280,65280,65280}
MAKE/FREE/O/D/N=10 green={32768,43520,54528,65280}
MAKE/FREE/O/D/N=10 blue={32768,32768,32768,32768}
j=0
wave/T elements
string selt=wave2list(elements)
string str2exec="",ListIn="",LocalListIn=""
for(k=0;k<nbbag;k+=1)
	if(SelBagOper[k]==48)
		sRatio="Ratio_"+ListBagOper[k]+"Oper"+num2str(quant)
		sMolbag="Molbag_"+ListBagOper[k]
		wave Molbag=$sMolbag
		nbmol=dimsize(Molbag,1)
		Make/O/D/N=(nbmol,2) $sRatio
		wave Ratio=$sRatio
		for(m=0;m<2;m+=1)
			str2exec=sRatio+"[]["+num2str(m)+"]="+parse2eval(ListRatioOper[m][1],ListBagOper[k])
			Execute str2exec
		endfor
		Appendtograph/W=$opername Ratio[][0] vs Ratio[][1]
		Modifygraph/W=$opername mode($sRatio)=3, marker($sRatio)=19, useMrkStrokeRGB($sRatio)=1, rgb($sRatio)=(red[j],green[j],blue[j])
		j+=1
	endif
endfor
Label/W=$opername bottom,ListRatioOper[1][1]
Label/W=$opername left,ListRatioOper[0][1]
end

function/S parse2eval(str,name)
string str,name
wave/T elements
string selt=wave2list(elements)
string LocalListIn=str
string res=""
variable V_flag=1,kpar
string elt1="",elt2="",elt3=""
string eltreplacor=""
for(kpar=0;strlen(LocalListIn) && V_flag;kpar+=1)
	eltreplacor=""
	Splitstring/E="([[:digit:]]*[.]*[[:digit:]]*[[:alpha:]]+)" LocalListIn, elt1
	Splitstring/E="([[:digit:]]+[.]*[[:digit:]]*)" elt1,elt2
	Splitstring/E="([[:alpha:]]+)" elt1,elt3
	if(strlen(elt2)==0)
		elt2="1"
	endif
	//print elt1,elt2,elt3, V_flag
	res+=LocalListIn[0,strsearch(LocalListIn,elt1,0)-1]
	LocalListIn=LocalListIn[strlen(elt1)+strsearch(LocalListIn,elt1,0),inf]
	if(V_flag)
        eltreplacor="("+elt2+"*"+"Molbag_"+name+"["+num2str(whichlistitem(elt3,selt,";",0,1))+"][x])"
        if(stringmatch(elt1,"mass"))
            eltreplacor=""+elt2+"*"+"Ciblebag_"+name+"[x][1]"
        elseif(stringmatch(elt1,"14N"))
            eltreplacor="Molbag_"+name+"[108][x]"
        elseif(stringmatch(elt1,"15N"))
            eltreplacor="Molbag_"+name+"[109][x]"
        elseif(stringmatch(elt1,"1H"))
            eltreplacor="Molbag_"+name+"[110][x]"
        elseif(stringmatch(elt1,"2H"))
            eltreplacor="Molbag_"+name+"[111][x]"
        elseif(stringmatch(elt1,"12C"))
            eltreplacor="Molbag_"+name+"[112][x]"
        elseif(stringmatch(elt1,"13C"))
            eltreplacor="Molbag_"+name+"[113][x]"
        endif
    else
		eltreplacor=LocalListIn
	endif
	res+=eltreplacor
endfor
return res
end



function InitKnariesLists(quant)
variable quant
string opername="AdvMolManOper"+num2str(quant)
if(wintype(opername)==0)
	return 0
endif
razTraces(opername)
string sKnaries, sMolbag
string sSelEltOper="SelEltOper"+num2str(quant)
string sListEltOper="ListEltOper"+num2str(quant)
string sSelBagOper="SelBagOper"+num2str(quant)
string sListBagOper="ListBagOper"+num2str(quant)
string sKnariesPoles="KnariesPolesOper"+num2str(quant)
string sIndxEltOper="IndxEltOper"+num2str(quant)
string sIndxTickedOper="IndxTickedOper"+num2str(quant)
string sPolesNames="PolesNamesOper"+num2str(quant)
wave/T ListBagOper=$sListBagOper
wave/T ListEltOper=$sListEltOper
wave IndxEltOper=$sIndxEltOper
wave SelBagOper=$sSelBagOper
wave SelEltOper=$sSelEltOper
Extract/O ListEltOper, $sPolesNames, SelEltOper==48
Extract/O IndxEltOper, $sIndxTickedOper, SelEltOper==48
wave KnariesPoles=$sKnariesPoles
wave IndxTickedOper=$sIndxTickedOper
wave/T PolesNames=$sPolesNames
variable dim=dimsize(PolesNames,0), nbbag=dimsize($sSelBagOper,0),nbmol,k,j
Duplicate/O placepoles(dim) $sKnariesPoles
Appendtograph/W=$opername KnariesPoles[][1] vs KnariesPoles[][0]
Appendtograph/W=$opername KnariesPoles[][1] vs KnariesPoles[][0]
Modifygraph/W=$opername mode($sKnariesPoles#0)=3,marker($sKnariesPoles#0)=5,msize($sKnariesPoles#0)=5, mirror=1,nticks=0,axOffset=-3,rgb($sKnariesPoles#0)=(0,0,0)
Modifygraph/W=$opername mode($sKnariesPoles#1)=3,marker($sKnariesPoles#1)=1, textMarker($sKnariesPoles#1)={$sPolesNames,"default",0,0,5,0.00,0.00}, rgb($sKnariesPoles#1)=(0,0,0)
MAKE/FREE/O/D/N=10 red={65280,65280,65280,65280}
MAKE/FREE/O/D/N=10 green={32768,43520,54528,65280}
MAKE/FREE/O/D/N=10 blue={32768,32768,32768,32768}
j=0
for(k=0;k<nbbag;k+=1)
	if(SelBagOper[k]==48)
	sKnaries="KnariesPos_"+ListBagOper[k]+"Oper"+num2str(quant)
	sMolbag="Molbag_"+ListBagOper[k]
	wave Molbag=$sMolbag
	nbmol=dimsize(Molbag,1)
	Make/O/FREE/D/N=(nbmol,dim) obs
	obs=Molbag[IndxTickedOper[y]][x]
	MatrixOP/FREE/O somme=sumrows(obs)
	obs/=somme[x]
	MatrixOP/FREE/O pos=obs x KnariesPoles
	Duplicate/O pos $sKnaries
	wave Knaries=$sKnaries
	Appendtograph/W=$opername Knaries[][1] vs Knaries[][0]
	Modifygraph/W=$opername mode($sKnaries)=3, marker($sKnaries)=19, useMrkStrokeRGB($sKnaries)=1, rgb($sKnaries)=(red[j],green[j],blue[j])
	j+=1
	endif
endfor
end

function/WAVE placepoles(n)
variable n
make/FREE/O/D/N=(n,2) poles=(y==x)/1
poles[][0]=cos(x*2*pi/n+pi/2)
poles[][1]=sin(x*2*pi/n+pi/2)
return poles
end


Function AdvMolManOperHook(s)
	STRUCT WMWinHookStruct &s
	Variable hookResult = 0
variable quant=str2num(replacestring("AdvMolManOper",s.winName,""))
string molsel="selmol"+num2str(quant)
wave selmol=$molsel
string panelname=replacestring("AdvMolManOper",s.winName,"advmolmanpanel")
string Opername=s.winName
variable hitpoint=str2num(stringbykey("HITPOINT", TraceFromPixel(s.mouseLoc.h, s.mouseLoc.v,""))),k
string tracename=stringbykey("TRACE", TraceFromPixel(s.mouseLoc.h, s.mouseLoc.v,""))
string infotrace=traceinfo(Opername,tracename,0)
string name=trace2bagname(tracename)
string tracelist=tracenamelist(Opername,";",1)
string hittracelist=""
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
for(k=0;k<itemsinlist(tracelist);k+=1)
	if(strsearch(stringfromlist(k,tracelist),name,0)!=-1)
		hittracelist+=stringfromlist(k,tracelist)+";"
	endif
endfor
if(strlen(name))
	//Tag/W=$Opername/C/N=text0/A=LB/Y=(5)/X=(5)/L=2 $tracename , hitpoint, listbag[hitpoint]+" "+num2str(ciblebag[hitpoint][3])
	MajTagAdvMolMan(quant,name,hitpoint)
endif
	switch(s.eventCode)
		case 2:
		CleanupOper(quant)
			break
		case 3:
		selmol[][0]=64
		selmol[hitpoint][0]=65
		listbox mollist win=$panelname, row=max(-1,hitpoint-2)
		if(s.eventMod==16 || s.eventMod==17)
			string scolor=replacestring("zColor(x)={",stringbykey("RECREATION",infotrace),"")
			scolor=replacestring(",*,*,Classification}",scolor,"")
			wave color=$scolor
			variable class=color[hitpoint]
			Extract/FREE/INDX/O color,ClassMemberIndx, color==class
			string item1="Copy this class to a new list;"
			//string item2="!!! Delete this class from "+name+" !!!;"
			popupcontextualmenu item1//+item2
				switch(V_flag)
					case 1:
						string newname
						prompt newname,"Copy this class to the following new list"
						Doprompt "Chose a destination list",newname
						if(1-V_flag)
							MolbagExtractor(name,ClassMemberIndx,newname)
							majbaglist(quant)
						endif
						break
					case 2:
						MolbagWithdrawFrom(ClassMemberIndx,name)
						majbaglist(quant)
						break
				endswitch
		endif
			break
		case 6:
		MajAdvMolManOperWindow(quant)
			break
		case 12:
			break
	endswitch
	return hookResult		// 0 if nothing done, else 1
end

function MajAdvMolManOperWindow(quant)
variable quant
string opername="AdvMolManOper"+num2str(quant)
getwindow $opername wsizeDC
variable width=v_right,height=v_bottom
controlinfo/W=$opername tab0
TabControl tab0, win=$opername, pos={1,1},size={width,v_height}
razcontrol(opername,"TabCon_*")
setcontrolsOPER(V_value,opername)
end

Function Pilottab(tca) : TabControl
	STRUCT WMTabControlAction &tca
	switch( tca.eventCode )
		case 2: // mouse up
			Variable tab = tca.tab
			razcontrol(tca.win,"TabCon_*")
			setcontrolsOPER(tab,tca.win)
			switch(tab)
				case 0:
					break
				case 1:
					break
			endswitch
			break
		case -1:
			break
	endswitch
	return 0
End

function razcontrol(win,str)//vire tous les controles qui matchent str dans la fenetre win
string win,str
string con2kill=ControlNameList(win,";",str),lecon
variable n=itemsinlist(con2kill),k
//print con2kill
for(k=0;k<n;k+=1)
	lecon=stringfromlist(k,con2kill)
	killcontrol/W=$win $lecon
endfor
end

function razTraces(graphname)
string graphname
string tracelist=tracenamelist(graphname,";",1)
string annotlist=annotationlist(graphname)
variable k,n=itemsinlist(tracelist)
string latrace,lannot
for(k=0;k<n;k+=1)
	latrace=stringfromlist(n-1-k,tracelist)
	removefromgraph/W=$graphname $latrace
endfor
n=itemsinlist(annotlist)
for(k=0;k<n;k+=1)
	lannot=stringfromlist(n-1-k,annotlist)
	textbox/W=$graphname/K/N=$lannot
endfor
end

function MajAdvMolManGraph(quant)
variable quant
string graphname="AdvMolManGraph"+num2str(quant)
if(wintype(graphname))
string panelname="AdvMolManPanel"+num2str(quant)
string molsel="selmol"+num2str(quant)
string mollist="mollist"+num2str(quant)
string bagsel="seladvmolman"+num2str(quant)
string baglist="advmolmanlist"+num2str(quant)
wave selWavebag=$bagsel
wave/T listWavebag=$baglist
wave/T elements
variable n=dimsize(selWavebag,0),k,asked,present,j,m,count,nbtr
string ongraph=tracenamelist(graphname,";",1), trace2kill="removefromgraph/W="+graphname+" "+tracenamelist(graphname,",",1)
trace2kill=trace2kill[0,strlen(trace2kill)-2]
string tracename=""
Execute trace2kill
ongraph=replacestring("Ciblebag_",ongraph,"")
ongraph=replacestring("Molbag_",ongraph,"")
string top="",bottom="",abscisse="",S_value=""
controlinfo/W=$graphname AbsMenu
abscisse=S_Value
controlinfo/W=$graphname TopMenu
top=S_Value
controlinfo/W=$graphname BottomMenu
bottom=S_Value
Make /O /N=14 rougeCompareAttrib={0,65000,0,0,33000,33000,60000,12500,0,0,46000,16000,33000,46000,65000}
Make /O /N=14 vertCompareAttrib={0,0,33000,16000,0,33000,49000,53000,0,65000,48000,33000,33000,3000,16000}
Make /O /N=14 bleuCompareAttrib={0,0,0,33000,33000,33000,5000,14000,65000,0,17000,33000,0,3000,0}
wave res=NonZeroElement(Ticked2List(selWavebag,listWavebag))
for(k=0;k<n;k+=1)
	asked=(selWavebag[k][5]==48 || mod(selWavebag[k][0],2))
	string lebag=listWavebag[k][1]
	string refcible="Ciblebag_"+lebag, refmolbag="Molbag_"+lebag
	tracename=""
	wave wavecible=$refcible
	wave wavemol=$refmolbag
	if(asked==1)
		string laxe="",testr=""
		tracename=""
		string legendstring=""
		for(j=0;j<2;j+=1)
			if(j==0)
				laxe="lefttop"
				testr=top
			else
				laxe="leftbottom"
				testr=bottom
			endif
			strswitch(testr)	// numeric switch
				case "Stoichiometry":
					m=dimsize(res,0)
					for(count=0;count<m;count+=1)
						AppendtoGraph/W=$graphname/L=$laxe wavemol[res[count]][] vs wavecible[][botmenu2idx(abscisse)]
						tracename=refmolbag+"#"+num2str(count)
						ModifyGraph/W=$graphname mode($tracename)=3, zmrkSize($tracename)={wavecible[*][2],*,*,1,10}, marker($tracename)=42
						ModifyGraph/W=$graphname rgb($tracename)=(rougeCompareAttrib[count],vertCompareAttrib[count],bleuCompareAttrib[count])
						legendstring+="\\s("+tracename+") "+elements[res[count]]+"\r"
					endfor
					legendstring=legendstring[0,strlen(legendstring)-2]
					Legend/W=$graphname/C/N=$laxe/J/X=-40.00/Y=-12.00/A=MC legendstring
					break
				default:
					count=0
					AppendtoGraph/W=$graphname/L=$laxe wavecible[][botmenu2idx(testr)] vs wavecible[][botmenu2idx(abscisse)]
					tracename=refcible+"#"+num2str(count)
					ModifyGraph/W=$graphname mode($tracename)=3, rgb($tracename)=(19968,19968,49920), marker($tracename)=0, msize($tracename)=0
					ModifyGraph/W=$graphname grid($laxe)=1,zero($laxe)=1
					Legend/W=$graphname/C/N=$laxe/K
					break
			endswitch
		endfor
		ModifyGraph/W=$graphname axisEnab(lefttop)={0.51,1},axisEnab(leftbottom)={0.01,0.49},axisEnab(bottom)={0.07,1}
	endif
endfor
endif
end

function botmenu2idx(str)
string str
	strswitch(str)	// string switch
		case "Target Mass":
			return 1
		case "Intensity":
			return 2
		case "Relative mass bias":
			return 3
		case "Y=f(X)":
			return 4
		case "Cost function":
			return 0
	endswitch
end

function/WAVE NonZeroElement(bagnamelist)//returns a n-wave contening the index of the n elements whose stoechoimetric coeff is never zero in the baglist
//bagnamelist must be ;-separated and contain the name of the molbag, not the references
string bagnamelist
variable n=itemsinlist(bagnamelist),k
string slebag=""
wave/T elements
make/O/D/N=(dimsize(elements,0)) res=0
for(k=0;k<n;k+=1)
	slebag="Molbag_"+stringfromlist(k,bagnamelist)
	wave lebag=$slebag
	matrixOP/O res=res+sumrows(abs(lebag))
endfor
res=x*res/res
wavetransform/O zapnans res
return res
end

function/S InitAdvMolManMenu(bagsel,baglist)
wave bagsel
wave/T baglist
wave leselm=NonZeroElement(Ticked2List(bagsel,baglist))
variable k, n=dimsize(leselm,0)
string res="Y=DBE;Y=RT;"
wave/T elements
	for(k=0;k<n;k+=1)
		res+="Y="+elements[leselm[k]]+";"
	endfor
res=res[0,strlen(res)-2]
return res
end

function/S Ticked2List(bagsel,baglist)//returns a string list of molbag in baglist ticked in the bagsel
wave bagsel
wave/T baglist
variable n=dimsize(bagsel,0),k
string res=""
for(k=0;k<n;k+=1)
	if(bagsel[k][5]==48 || mod(bagsel[k][0],2))
		res+=baglist[k][1]+";"
	endif
endfor
return res
end

Function Wessimple(w,M) : FitFunc
	Wave w
	Variable M

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(M) = K/M * 1/(sqrt(2*pi*L*(1-P))) * exp( -0.5* (ln(M/b)/ln(a)-L)^2/ (L*(1-p)) )
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ M
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = L
	//CurveFitDialog/ w[1] = P
	//CurveFitDialog/ w[2] = K
	//CurveFitDialog/ w[3] = a
	//CurveFitDialog/ w[4] = b

	return w[2]/M * 1/(sqrt(2*pi*w[0]*(1-w[1]))) * exp( -0.5* (ln(M/w[4])/ln(w[3])-w[0])^2/ (w[0]*(1-w[1])) )
End

function unsat(mol)
wave mol
return 1+mol[5]+mol[6]/2-mol[0]/2
end

function/wave condirez()
wave W_coef,W_sigma, roi1,roi0
nvar V_chisq, V_npnts
duplicate/O roi1 transavg
Duplicate/O/FREE roi1 entropy, residue
residue=(roi1[x]-Wesslau(W_coef,roi0[x]))^2
variable somme=sum(roi1)
entropy=roi1/somme
entropy=entropy[x]*ln(entropy[x])
transavg=roi1*roi0
make/O/N=(numpnts(W_coef)+numpnts(W_coef)+4) condi_fitrez
condi_fitrez[0,4]=W_coef
condi_fitrez[5,9]=W_sigma[x-5]
condi_fitrez[10]=sqrt(sum(residue))/somme
condi_fitrez[11]=numpnts(roi0)
condi_fitrez[12]=sum(transavg)/sum(roi1)//massmoyenne
condi_fitrez[13]=sum(entropy)
condi_fitrez[0]=somme
killwaves transavg
return condi_fitrez
end

Window composite_article() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(316.5,42.5,1120.75,740) segzdm vs segzm
	AppendToGraph simu_probadm vs simu_mass
	AppendToGraph roiy vs roix
	AppendToGraph yplode vs xplode
	AppendToGraph cauchyyplode vs cauchyxplode
	AppendToGraph defmasstree vs tree
	AppendToGraph wave2stripped vs wave0stripped
	AppendToGraph incoxdm vs incox
	AppendToGraph roi2 vs roi0
	AppendToGraph currentseg[*][0] vs currentseg[*][1]
	AppendToGraph/L=Intensity wave1stripped vs wave0stripped
	AppendToGraph/L=Intensity roi1 vs roi0
	ModifyGraph margin(right)=113,gbRGB=(56576,56576,56576)
	ModifyGraph mode(simu_probadm)=3,mode(defmasstree)=4,mode(wave2stripped)=3,mode(incoxdm)=3
	ModifyGraph mode(roi2)=4,mode(wave1stripped)=1,mode(roi1)=4
	ModifyGraph marker(simu_probadm)=18,marker(defmasstree)=55,marker(wave2stripped)=19
	ModifyGraph marker(incoxdm)=42,marker(roi2)=19,marker(wave1stripped)=9,marker(roi1)=19
	ModifyGraph lSize(yplode)=1.2,lSize(cauchyyplode)=1.2,lSize(roi2)=2,lSize(currentseg)=2
	ModifyGraph lSize(wave1stripped)=1.5,lSize(roi1)=2
	ModifyGraph rgb(segzdm)=(65535,65535,65535),rgb(simu_probadm)=(0,52224,52224),rgb(yplode)=(0,0,65280)
	ModifyGraph rgb(defmasstree)=(34816,34816,34816),rgb(wave2stripped)=(0,0,0),rgb(incoxdm)=(26368,0,52224)
	ModifyGraph rgb(roi2)=(0,0,65280),rgb(currentseg)=(0,52224,52224),rgb(roi1)=(0,0,65280)
	ModifyGraph msize(defmasstree)=7,msize(wave2stripped)=1.1,msize(incoxdm)=5,msize(wave1stripped)=1.5
	ModifyGraph mrkThick(defmasstree)=1,mrkThick(wave2stripped)=0.01,mrkThick(roi2)=1
	ModifyGraph mrkThick(roi1)=1
	ModifyGraph opaque(incoxdm)=1,opaque(roi2)=1,opaque(roi1)=1
	ModifyGraph hideTrace(segzdm)=1,hideTrace(simu_probadm)=1,hideTrace(roiy)=1,hideTrace(yplode)=1
	ModifyGraph hideTrace(cauchyyplode)=1,hideTrace(defmasstree)=1,hideTrace(incoxdm)=1
	ModifyGraph hideTrace(currentseg)=1
	ModifyGraph useMrkStrokeRGB(simu_probadm)=1,useMrkStrokeRGB(defmasstree)=1,useMrkStrokeRGB(roi2)=1
	ModifyGraph useMrkStrokeRGB(roi1)=1
	ModifyGraph mrkStrokeRGB(simu_probadm)=(65535,65535,65535),mrkStrokeRGB(defmasstree)=(65535,65535,65535)
	ModifyGraph mrkStrokeRGB(wave2stripped)=(65535,65535,65535),mrkStrokeRGB(roi2)=(65535,65535,65535)
	ModifyGraph mrkStrokeRGB(roi1)=(65535,65535,65535)
	ModifyGraph zmrkSize(defmasstree)={treeerror,*,10,6,0}
	ModifyGraph zColor(defmasstree)={treeerror,0,10,Bathymetry9,1},zColor(wave2stripped)={logwave1,*,*,YellowHot}
	ModifyGraph zColor(wave1stripped)={wave1stripped,*,*,YellowHot}
	ModifyGraph zColorMin(wave2stripped)=(65535,65535,65535)
	ModifyGraph logZColor(wave1stripped)=1
	ModifyGraph grid(bottom)=1
	ModifyGraph log(Intensity)=1
	ModifyGraph zero(left)=1
	ModifyGraph mirror=1
	ModifyGraph nticks(left)=10,nticks(bottom)=40
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(left)=16
	ModifyGraph standoff(bottom)=0
	ModifyGraph axOffset(bottom)=-0.238095
	ModifyGraph axThick=2
	ModifyGraph gridRGB(bottom)=(65535,65535,65535)
	ModifyGraph gridStyle(bottom)=2
	ModifyGraph lblPos(left)=69,lblPos(Intensity)=66
	ModifyGraph lblLatPos(left)=-6,lblLatPos(Intensity)=-4
	ModifyGraph freePos(Intensity)=0
	ModifyGraph axisEnab(left)={0,0.49}
	ModifyGraph axisEnab(Intensity)={0.51,1}
	Label left "Mass defect"
	Label bottom "m/z"
	Label Intensity "Intensity (log scale)"
	SetAxis left -0.5,0.5
	SetAxis bottom 150,650
	SetAxis Intensity*,5701076.4
	ColorScale/C/N=text0/F=0/A=MC/X=57.36/Y=25.40 trace=wave2stripped, heightPct=49
	ColorScale/C/N=text0 nticks=10
	AppendText "\\f01Colorscale (log\\B10\\M units)"
	Legend/C/N=text1/J/A=MC/X=24.19/Y=41.32 "\\Z24A\r\\s(roi1) roi1"
	Legend/C/N=text1_1/J/A=MC/X=-11.75/Y=-35.93 "\\Z24B\r\\s(roi1) roi1"
	ShowTools/A
	SetWindow kwTopWin,hook=getdmvmMSMSmode,hookevents=3
EndMacro

function histomurch()
wave massmurch, roi0,nbmurch
variable n=numpnts(roi0),k
for(k=0;k<n;k+=1)
	nbmurch[floor(roi0[k])-150]+=1
endfor
end

Window dmvm_1() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(391.5,110,1172.25,657.5) simu_probadm vs simu_mass
	AppendToGraph roiy vs roix
	AppendToGraph yplode vs xplode
	AppendToGraph cauchyyplode vs cauchyxplode
	AppendToGraph defmasstree vs tree
	AppendToGraph incoxdm vs incox
	AppendToGraph roi2 vs roi0
	AppendToGraph currentseg[*][0] vs currentseg[*][1]
	AppendToGraph simu_alkanes_2stripped vs simu_alkanes_0stripped
	AppendToGraph DefmassForm_1 vs Formule_1mass
	AppendToGraph DefmassForm_2 vs Formule_2mass
	AppendToGraph DefmassForm_3 vs Formule_3mass
	AppendToGraph DefmassForm_4 vs Formule_4mass
	AppendToGraph DefmassForm_5 vs Formule_5mass
	AppendToGraph DefmassForm_6 vs Formule_6mass
	AppendToGraph DefmassForm_7 vs Formule_7mass
	AppendToGraph DefmassForm_8 vs Formule_8mass
	AppendToGraph DefmassForm_9 vs Formule_9mass
	AppendToGraph DefmassForm_10 vs Formule_10mass
	AppendToGraph DefmassForm_11 vs Formule_11mass
	AppendToGraph CH2streak4_2stripped vs CH2streak4_0stripped
	AppendToGraph CH2streak5_2stripped vs CH2streak5_0stripped
	AppendToGraph CH2streak6_2stripped vs CH2streak6_0stripped
	AppendToGraph wave2stripped vs wave0stripped
	ModifyGraph margin(top)=45,margin(right)=170,gbRGB=(47872,47872,47872)
	ModifyGraph mode(simu_probadm)=3,mode(defmasstree)=3,mode(incoxdm)=3,mode(roi2)=3
	ModifyGraph mode(simu_alkanes_2stripped)=4,mode(DefmassForm_1)=3,mode(DefmassForm_2)=3
	ModifyGraph mode(DefmassForm_3)=3,mode(DefmassForm_4)=3,mode(DefmassForm_5)=3,mode(DefmassForm_6)=3
	ModifyGraph mode(DefmassForm_7)=3,mode(DefmassForm_8)=3,mode(DefmassForm_9)=3,mode(DefmassForm_10)=3
	ModifyGraph mode(DefmassForm_11)=3,mode(CH2streak4_2stripped)=4,mode(CH2streak5_2stripped)=4
	ModifyGraph mode(CH2streak6_2stripped)=4,mode(wave2stripped)=3
	ModifyGraph marker(simu_probadm)=18,marker(defmasstree)=55,marker(incoxdm)=42,marker(roi2)=19
	ModifyGraph marker(simu_alkanes_2stripped)=19,marker(DefmassForm_1)=18,marker(DefmassForm_2)=18
	ModifyGraph marker(DefmassForm_3)=18,marker(DefmassForm_4)=18,marker(DefmassForm_5)=18
	ModifyGraph marker(DefmassForm_6)=18,marker(DefmassForm_7)=18,marker(DefmassForm_8)=18
	ModifyGraph marker(DefmassForm_9)=18,marker(DefmassForm_10)=18,marker(DefmassForm_11)=18
	ModifyGraph marker(CH2streak4_2stripped)=13,marker(CH2streak5_2stripped)=41,marker(CH2streak6_2stripped)=40
	ModifyGraph marker(wave2stripped)=19
	ModifyGraph lSize(yplode)=1.2,lSize(cauchyyplode)=1.2,lSize(roi2)=1.2,lSize(currentseg)=2
	ModifyGraph lSize(simu_alkanes_2stripped)=2,lSize(CH2streak4_2stripped)=2,lSize(CH2streak5_2stripped)=2
	ModifyGraph lSize(CH2streak6_2stripped)=2
	ModifyGraph lStyle(simu_alkanes_2stripped)=3,lStyle(CH2streak4_2stripped)=3,lStyle(CH2streak5_2stripped)=3
	ModifyGraph lStyle(CH2streak6_2stripped)=3
	ModifyGraph rgb(simu_probadm)=(0,52224,52224),rgb(yplode)=(0,0,65280),rgb(defmasstree)=(34816,34816,34816)
	ModifyGraph rgb(incoxdm)=(26368,0,52224),rgb(roi2)=(36864,14592,58880),rgb(currentseg)=(0,52224,52224)
	ModifyGraph rgb(simu_alkanes_2stripped)=(0,52224,52224),rgb(DefmassForm_1)=(0,130,1000)
	ModifyGraph rgb(DefmassForm_2)=(0,2323,1000),rgb(DefmassForm_3)=(0,208,1000),rgb(DefmassForm_4)=(0,1591,1000)
	ModifyGraph rgb(DefmassForm_5)=(0,1047,1000),rgb(DefmassForm_6)=(0,1909,1000),rgb(DefmassForm_7)=(0,4419,1000)
	ModifyGraph rgb(DefmassForm_8)=(0,26,1000),rgb(DefmassForm_9)=(0,652,1000),rgb(DefmassForm_10)=(0,339,1000)
	ModifyGraph rgb(DefmassForm_11)=(0,233,1000),rgb(CH2streak4_2stripped)=(0,0,0),rgb(CH2streak5_2stripped)=(0,0,0)
	ModifyGraph rgb(CH2streak6_2stripped)=(0,0,0),rgb(wave2stripped)=(0,0,0)
	ModifyGraph msize(defmasstree)=7,msize(incoxdm)=5,msize(roi2)=1.5,msize(simu_alkanes_2stripped)=4
	ModifyGraph msize(CH2streak4_2stripped)=8,msize(CH2streak5_2stripped)=8,msize(CH2streak6_2stripped)=8
	ModifyGraph msize(wave2stripped)=5
	ModifyGraph mrkThick(defmasstree)=1,mrkThick(simu_alkanes_2stripped)=2,mrkThick(CH2streak4_2stripped)=3
	ModifyGraph mrkThick(CH2streak5_2stripped)=3,mrkThick(CH2streak6_2stripped)=3,mrkThick(wave2stripped)=0.01
	ModifyGraph opaque(incoxdm)=1
	ModifyGraph hideTrace(simu_probadm)=1,hideTrace(yplode)=1,hideTrace(cauchyyplode)=1
	ModifyGraph useMrkStrokeRGB(simu_probadm)=1,useMrkStrokeRGB(defmasstree)=1,useMrkStrokeRGB(roi2)=1
	ModifyGraph useMrkStrokeRGB(simu_alkanes_2stripped)=1,useMrkStrokeRGB(DefmassForm_1)=1
	ModifyGraph useMrkStrokeRGB(DefmassForm_2)=1,useMrkStrokeRGB(DefmassForm_3)=1,useMrkStrokeRGB(DefmassForm_4)=1
	ModifyGraph useMrkStrokeRGB(DefmassForm_5)=1,useMrkStrokeRGB(DefmassForm_6)=1,useMrkStrokeRGB(DefmassForm_7)=1
	ModifyGraph useMrkStrokeRGB(DefmassForm_8)=1,useMrkStrokeRGB(DefmassForm_9)=1,useMrkStrokeRGB(DefmassForm_10)=1
	ModifyGraph useMrkStrokeRGB(DefmassForm_11)=1,useMrkStrokeRGB(CH2streak4_2stripped)=1
	ModifyGraph useMrkStrokeRGB(CH2streak5_2stripped)=1,useMrkStrokeRGB(CH2streak6_2stripped)=1
	ModifyGraph mrkStrokeRGB(simu_probadm)=(65535,65535,65535),mrkStrokeRGB(wave2stripped)=(65535,65535,65535)
	ModifyGraph zmrkSize(defmasstree)={treeerror,*,10,6,0},zmrkSize(wave2stripped)={logwave1,*,*,1,5}
	ModifyGraph zColor(defmasstree)={treeColor,*,*,Web216,1},zColor(wave2stripped)={logwave1,4.31215111173074,*,YellowHot}
	ModifyGraph zColorMin(wave2stripped)=(65535,65535,65535)
	ModifyGraph zero(left)=1
	ModifyGraph mirror=1
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	Label left "\\Z16Mass defect"
	Label bottom "\\Z16Mass"
	SetAxis left 0.235787841728891,0.409687486239092
	SetAxis bottom 300.36519,398
	ShowInfo
	Tag/C/N=text01/D=1.2/X=-13.03/Y=14.14/L=1/TL={lThick=1.2} simu_alkanes_2stripped, 22, "\\JCAlkanes guideline:\rH\\B2\\M \\f02\\K(0,0,65280)(CH\\B2\\M)\\Bn"
	Tag/C/N=text02/D=1.2/A=LC/X=16.00/Y=-8.00/L=1/TL={lThick=1.2} CH2streak4_2stripped, 17, "H\\S+\\M + (NH)C\\B\\f04\\f053\\M\\f02 \\K(0,0,65280)(CH\\B2\\M)\\Bn"
	Tag/C/N=text03/D=1.2/A=LC/X=18.00/Y=-9.00/L=1/TL={lThick=1.2} CH2streak5_2stripped, 16, "H\\S+\\M + (NH)C\\B\\f04\\f054\\M\\f02 \\K(0,0,65280)(CH\\B2\\M)\\Bn"
	Tag/C/N=text04/D=1.2/A=LC/X=20.00/Y=-10.00/L=1/TL={lThick=1.2} CH2streak6_2stripped, 15, "H\\S+\\M + (NH)C\\B\\f04\\f055\\M\\f02 \\K(0,0,65280)(CH\\B2\\M)\\Bn"
	ValDisplay valdisp0,pos={6,1},size={89,14},title="ROI size"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"numpnts(roi1)"
	ValDisplay valdisp1,pos={106,1},size={99,14},title="IntCum %"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"round(100*(sum(roi1)/sum(wave1stripped)))"
	ValDisplay valdisp2,pos={326,1},size={74,14},title="InPnts"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000},value= #"twosig"
	SetVariable setvar0,pos={205,1},size={120,16},proc=Proc_setintenscumthres,title="IntCumThres"
	SetVariable setvar0,limits={0,1,0.1},value= intenscumthres
	SetVariable setvar1,pos={6,18},size={113,16},proc=SetVarProc_6,title="Depth of tree"
	SetVariable setvar1,limits={0,inf,1},value= _NUM:1
	ValDisplay valdisp3,pos={124,18},size={142,14},title="Segment length :"
	ValDisplay valdisp3,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp3,value= #"sqrt((currentseg[0][1]-currentseg[1][1])^2+(currentseg[0][0]-currentseg[1][0])^2)"
	PopupMenu popup0,pos={269,18},size={123,21},proc=PopMenuProc_3,title="Use this list:"
	PopupMenu popup0,mode=2,popvalue="Factors",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	SetWindow kwTopWin,hook=getdmvmMSMSmode,hookevents=3
EndMacro

function killmissfit(seuil)
variable seuil
variable k, n=itemsinlist(wavelist("molecule_*",";",""))
string lamol, lafor, blaz, ledef, lamas
wave/Z trans, index_formules,list_targets
wave/T list_formules
duplicate/O index_formules color statmag statsou statox statsod
make/O/D/N=(n) transdelta
transdelta=list_targets[x][3]
	for(k=n;k>0;k-=1)
		lamol="molecule_"+num2str(k)
		lamas="formule_"+num2str(k)+"mass"
		lafor="formule_"+num2str(k)
		ledef="DefmassForm_"+num2str(k)
		duplicate/O $lamol trans
		trans=abs(trans)
		if(abs(transdelta[k-1])>seuil)
		Removefromgraph/W=agregateur $lafor
		Removefromgraph/W=dmvm $ledef
		Killwaves/Z $lafor,$lamas, $lamol, $ledef
		deletepoints/M=0 k-1, 1, list_formules, index_formules, list_targets
		MolbagWithdrawFrom({k-1},"Current")
		endif
	endfor
	majformulewavesnames()
	majliste()
	majlegende()
	killwaves transdelta
end
macro tempo_rose()
	appendtograph rosepentesYselect vs rosepentesXselect
	AppendToGraph roseDesPentesFacteur0Y vs roseDesPentesFacteur0X
	AppendToGraph roseDesPentesFacteur1Y vs roseDesPentesFacteur1X
	AppendToGraph roseDesPentesFacteur2Y vs roseDesPentesFacteur2X
	AppendToGraph roseDesPentesFacteur3Y vs roseDesPentesFacteur3X
	AppendToGraph roseDesPentesFacteur4Y vs roseDesPentesFacteur4X
	AppendToGraph roseDesPentesFacteur5Y vs roseDesPentesFacteur5X
	AppendToGraph roseDesPentesFacteur6Y vs roseDesPentesFacteur6X
	AppendToGraph roseDesPentesFacteur7Y vs roseDesPentesFacteur7X
	AppendToGraph roseDesPentesFacteur8Y vs roseDesPentesFacteur8X
	AppendToGraph roseDesPentesFacteur9Y vs roseDesPentesFacteur9X
	AppendToGraph roseDesPentesFacteur10Y vs roseDesPentesFacteur10X
	AppendToGraph roseDesPentesFacteur11Y vs roseDesPentesFacteur11X
	ModifyGraph lSize(rosePentesYSelect)=2

	ModifyGraph rgb(roseDesPentesFacteur1Y)=(0,0,0),rgb(roseDesPentesFacteur2Y)=(0,0,0)
	ModifyGraph rgb(roseDesPentesFacteur3Y)=(0,0,0),rgb(roseDesPentesFacteur4Y)=(0,0,0)
	ModifyGraph rgb(roseDesPentesFacteur5Y)=(0,0,0),rgb(roseDesPentesFacteur6Y)=(0,0,0)
	ModifyGraph rgb(roseDesPentesFacteur7Y)=(0,0,0),rgb(roseDesPentesFacteur8Y)=(0,0,0)
	ModifyGraph rgb(roseDesPentesFacteur9Y)=(0,0,0),rgb(roseDesPentesFacteur10Y)=(0,0,0)
	ModifyGraph rgb(roseDesPentesFacteur11Y)=(0,0,0)

	ModifyGraph muloffset(roseDesPentesFacteur0Y)={0.5,0.5},muloffset(roseDesPentesFacteur1Y)={0.5,0.5}
	ModifyGraph muloffset(roseDesPentesFacteur2Y)={0.5,0.5},muloffset(roseDesPentesFacteur3Y)={0.5,0.5}
	ModifyGraph muloffset(roseDesPentesFacteur4Y)={0.5,0.5},muloffset(roseDesPentesFacteur5Y)={0.5,0.5}
	ModifyGraph muloffset(roseDesPentesFacteur6Y)={0.5,0.5},muloffset(roseDesPentesFacteur7Y)={0.5,0.5}
	ModifyGraph muloffset(roseDesPentesFacteur8Y)={0.5,0.5},muloffset(roseDesPentesFacteur9Y)={0.5,0.5}
	ModifyGraph muloffset(roseDesPentesFacteur10Y)={0.5,0.5},muloffset(roseDesPentesFacteur11Y)={0.5,0.5}
	Tag/C/N=roseDesPentesFacteur0Y roseDesPentesFacteur0Y, 1, "H\\B2\\MC\\B1\\M"
	Tag/C/N=roseDesPentesFacteur1Y roseDesPentesFacteur1Y, 1, "H\\B8\\MC\\B5\\M"
	Tag/C/N=roseDesPentesFacteur2Y roseDesPentesFacteur2Y, 1, "H\\B2\\MC\\B2\\M"
	Tag/C/N=roseDesPentesFacteur3Y/X=9.36/Y=3.46 roseDesPentesFacteur3Y, 1, "H\\B6\\MC\\B5\\M"
	Tag/C/N=roseDesPentesFacteur4Y/X=-3.47/Y=6.29 roseDesPentesFacteur4Y, 1, "H\\B6\\MC\\B2\\M"
	Tag/C/N=roseDesPentesFacteur5Y/X=1.73/Y=9.59 roseDesPentesFacteur5Y, 1, "H\\B8\\MC\\B3\\M"
	Tag/C/N=roseDesPentesFacteur6Y roseDesPentesFacteur6Y, 1, "H\\B12\\MC\\B5\\M"
	Tag/C/N=roseDesPentesFacteur7Y/X=12.48/Y=2.20 roseDesPentesFacteur7Y, 1, "C\\B1\\M"
	Tag/C/N=roseDesPentesFacteur8Y/X=8.32/Y=3.30 roseDesPentesFacteur8Y, 1, "H\\B4\\MC\\B3\\M"
	Tag/C/N=roseDesPentesFacteur9Y/X=12.48/Y=3.30 roseDesPentesFacteur9Y, 1, "H\\B2\\MC\\B3\\M"
	Tag/C/N=roseDesPentesFacteur10Y/X=13.34/Y=2.99 roseDesPentesFacteur10Y, 1, "H\\B2\\MC\\B2\\MO\\B1\\M"
	Tag/C/N=roseDesPentesFacteur11Y/X=1.04/Y=10.38 roseDesPentesFacteur11Y, 1, "H\\B2\\M"
EndMacro

Structure Molbag
	string name
	WAVE Target
	WAVE/T List
	WAVE Matrix
	WAVE Charge
	variable size
EndStructure

Function MakeItBag(name)
string name
STRUCT Molbag lebag
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
if(waveexists($smolbag) && waveexists($schargebag) && waveexists($slistbag) && waveexists($sciblebag))
	lebag.Matrix=molbag
	lebag.Charge=chargebag
	lebag.List=listbag
	lebag.Target=ciblebag
	lebag.size=dimsize(lebag.Matrix,1)
endif
return lebag.size
end

Function belongto(mol,charge,nom,cible,bag,criter)//N pour la formule, E<XXX pour l'erreur inférieure à XXX
wave mol,cible
variable charge
string nom, bag,criter
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable res=0
variable k=0,n=dimsize(molbag,1), j=itemsinlist(criter),m
	for(k=0;k<n && res==0;k+=1)
		res=1
		Make/O/FREE/D/N=(114) mol2 = molbag[x][k]
		variable charge2=chargebag[k]
		string nom2=listbag[k]
		Make/O/FREE/D/N=5 cible2 = ciblebag[k][x]
		for(m=0;m<j;m+=1)
		string lecrit=stringfromlist(m,criter)
			strswitch(lecrit)
				case "N":
					res*=stringmatch(nom,nom2)
					break
				case "E<*":
					res*=abs(cible2[3])<str2num(lecrit[2,strlen(lecrit)-1])
					break
				default:
					res=0
			endswitch
		endfor
	endfor
print mol,charge,nom,cible,k,m,str2num(lecrit[2,strlen(lecrit)-1])
return res
end

function/WAVE BagErrorMatchIndex(bag,thresh)//Returns the index of which molecule matches
string bag
variable thresh
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1)
make/FREE/O/D/N=(n) lescibles=abs(ciblebag[x][3])
extract/FREE/indx/O lescibles, index, lescibles<thresh
return index
end

function MolbagExtractor(bag,index,newbagname)// extract the indexth molecules to the newbagname NON DESIGNED
string bag, newbagname
wave index
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=numpnts(index),k
for(k=0;k<n;k+=1)
	Make/O/FREE/D/N=(114) mol2 = molbag[x][index[k]]
	variable charge2=chargebag[index[k]]
	string nom2=listbag[index[k]]
	Make/O/FREE/D/N=5 cible2 = ciblebag[index[k]][x]
	MolbagAddTo(mol2,charge2,nom2,cible2,newbagname)
endfor
end

function/wave Molbag2FullIndx(bagname)
string bagname
string schargebag="Chargebag_"+bagname
wave chargebag=$schargebag
duplicate/FREE chargebag res
res=x
return res
end

function Molbag2Roi(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1)
Make/O/D/N=(n) roi1,roi0,roi2
roi0=ciblebag[x][1]/(1-ciblebag[x][3]/1e6)
roi1=ciblebag[x][2]
roi2=defmass(roi0)
end

function/Wave Molbag2Table(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1)
Make/O/T/N=1 titles={"Measured Mass (u)","Exact Mass (u)","Relative Mass Bias (ppm u)","Intensity","Y_val"}
Duplicate/FREE/O NonZeroElement(bag) eltz
variable m2=numpnts(eltz), m1=numpnts(titles)
Make/FREE/O/T/N=(m2) eltzname
wave/T elements
insertpoints m1,m2, titles
titles[m1,*]=elements[eltz[x-m1]]
Make/FREE/O/N=(n,m1+m2) res
res[][0]=(ciblebag[x][1])
res[][1]=(ciblebag[x][1]/(1-1e-6*ciblebag[x][3]))
res[][2]=(ciblebag[x][3])
res[][3]=(ciblebag[x][2])
res[][4]=(ciblebag[x][4])
res[][m1,*]=(Molbag[eltz[y-m1]][x])
return res
end

function interfero(name1,name2)
string name1,name2
string sint1, sint2, smass1, smass2
sint1=name1+"_1stripped"
sint2=name2+"_1stripped"
smass1=name1+"_0stripped"
smass2=name2+"_0stripped"
wave int1=$sint1, int2=$sint2, mass1=$smass1, mass2=$smass2
variable n=numpnts(mass1), m=numpnts(mass2),k,j
make/O/D/N=(n,m) comp
comp=abs(mass1[x]-mass2[y])
Make/O/D/N=0 roi0,roi1,roi2
variable sum1=sum(int1), sum2=sum(int2)
for(k=0;k<n;k+=1)
	for(j=0;j<m;j+=1)
		if(comp[k][j]<1e-5)
			insertpoints 0,1, roi0,roi1,roi2
			roi0[0]=mass1[k]
			roi1[0]=(int1[k]/sum1-int2[j]/sum2)//(int1[k]/sum1+int2[j]/sum2)
			roi2[0]=defmass(roi0[0])
			j=m
		endif
	endfor
endfor
end

function/WAVE MolbagAvg(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1),m=dimsize(molbag,0),k, somme=0
Make/O/D/N=(m) res=0
for(k=0;k<n;k+=1)
	res+=ciblebag[k][2]*molbag[x][k]
	somme+=ciblebag[k][1]
endfor
res/=somme
return res
end

Window Graph3_2() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(950.25,280.25,1318.5,654.5) Formule_1 vs Formule_1mass
	AppendToGraph Formule_2 vs Formule_2mass
	ModifyGraph wbRGB=(34816,34816,34816),gbRGB=(47872,47872,47872)
	ModifyGraph mode=8
	ModifyGraph marker=19
	ModifyGraph rgb(Formule_1)=(0,0,65280),rgb(Formule_2)=(65280,0,0)
	ModifyGraph useMrkStrokeRGB=1
	ModifyGraph log(left)=1
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	ModifyGraph axRGB=(65535,65535,65535)
	ModifyGraph tlblRGB=(65535,65535,65535)
	ModifyGraph alblRGB=(65535,65535,65535)
	Label left "Normalized probability"
	Label bottom "Mass"
	SetAxis left 0.001,1
	SetAxis bottom 799.61998,804.50161
EndMacro

Window Graph6() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(722.25,76.25,1116.75,284.75) ratio1312 vs meanErr
EndMacro

function spanIsotope(numato,nbadd,listmasse)//aout2017
variable numato,nbadd
wave listmasse
wave mendev_masses, wave1stripped,wave0stripped
variable n=dimsize(listmasse,0),k
Make/O/D/N=(n,4) res_iso
Make/O/D/N=4 uplet, mass2check, error
for(k=0;k<n;k+=1)
	//pour chaque masse à verifier on construit le 4-uplet de masses à verifier
	mass2check={listmasse[k],listmasse[k]+mendev_masses[numato][1]-mendev_masses[numato][0],listmasse[k]+nbadd*mendev_masses[numato][0],listmasse[k]+nbadd*mendev_masses[numato][0]+mendev_masses[numato][1]-mendev_masses[numato][0]}
	uplet=wave1stripped[prochepic(mass2check,0)]
	error=abs(1e6*(mass2check-wave0stripped[prochepic(mass2check,0)])/mass2check)
	res_iso[k][0]=(uplet[3]/uplet[2]-uplet[1]/uplet[0])/nbadd
	res_iso[k][1]=listmasse[k]
	res_iso[k][2]=mean(error)
	res_iso[k][3]=mean(uplet)
endfor
end

function startGraph(inc,ref)//septembre2017
variable inc,ref
wave slopes, intercepts,diffm,dataindex,modules,currentseg
duplicate/O diffm costfunction
variable n=dimsize(costfunction,0),pente,k,j,curr,score,cont,test=0
pente=defmass(ref)/ref
costfunction=sqrt((slopes-pente)^2+(intercepts-(defmass(inc)-inc*pente))^2+(abs(diffm)-ref)^2)
make/O/D/N=1 roipnts=prochepic(inc,0), roiscore=0
for(j=0;test==0;j+=1)
	findvalue/V=(roipnts[numpnts(roipnts)-1]) dataindex
	curr=V_value
	score=inf
	cont=nan
	for(k=0;k<n;k+=1)
		findvalue/V=(costfunction[k][curr]) roiscore
		if(costfunction[k][curr]<score && v_value==-1)
			score=costfunction[k][curr]
			cont=k
		endif
	endfor
	test=(abs(mean(roiscore)-score)/sqrt(variance(roiscore))>30)
	//print abs(mean(roiscore)-score), abs(mean(roiscore)-score)/sqrt(variance(roiscore))
	insertpoints numpnts(roipnts),1, roipnts,roiscore
	roipnts[numpnts(roipnts)-1]=dataindex[cont]
	roiscore[numpnts(roipnts)-1]=score
	roipnts2roiz()
	doupdate
endfor
end

function pnt2index(n)
variable n
wave dataindex
findvalue/V=(n) dataindex
return V_value
end

function courbplan(pos,lesx,lesy)//renvoie la courbure plane des trois points de pos
wave pos,lesx,lesy
if(numpnts(pos)==3)
	variable xp,xs,yp,ys,courbure
	xp=(lesx[pos[2]]-lesx[pos[0]])/2
	yp=(lesy[pos[2]]-lesy[pos[0]])/2
	xs=lesx[pos[2]]-2*lesx[pos[1]]+lesx[pos[0]]
	ys=lesy[pos[2]]-2*lesy[pos[1]]+lesy[pos[0]]
	courbure=(xp*ys-xs*yp)/((xp^2+yp^2)^(3/2))
	return courbure
endif
end

Function MolbagDBE2fifth(name)//janvier2018
string name
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
ciblebag[][4]=1-molbag[0][x]/2+molbag[5][x]+molbag[6][x]/2
end

Function MolbagElement2fifth(name,element)
string name
variable element
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
ciblebag[][4]=molbag[element][x]
end

function MolbagResort(name,drive)
string name
variable drive
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
Duplicate/O chargebag index, wave2index
variable size=dimsize(index,0)
	switch(drive)	// numeric switch
		case 3: //Target
			wave2index=ciblebag[x][1]
			makeindex wave2index, index
			break
		case 4://Bias
			wave2index=abs(ciblebag[x][3])
			makeindex wave2index, index
			break
		case 5://intensity
			wave2index=ciblebag[x][2]
			makeindex wave2index, index
			break
		case 7: //Fifth col
			wave2index=ciblebag[x][4]
			makeindex wave2index, index
			break
		default:
			Killwaves wave2index, index
			return 0
			break
	endswitch
	if(index[0]<index[size-1])
		wavetransform/O flip index
	endif
Duplicate/O Chargebag newcharge
Duplicate/O Ciblebag newcible
Duplicate/O/T Listbag newlist
Duplicate/O Molbag newmol
multithread newmol=molbag[x][index[y]]
newlist=listbag[index[x]]
multithread newcharge=chargebag[index[x]]
multithread newcible=ciblebag[index[x]][y]
Duplicate/O newcharge Chargebag
Duplicate/O newcible Ciblebag
Duplicate/O/T newlist Listbag
Duplicate/O newmol Molbag
Killwaves newcharge,newcible,newlist,newmol,wave2index, index
end

Function ClassBiasPoly(name,Order,nbClust,quant)
string name
variable Order,nbClust,quant
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
string sClass="Clu_Cla9_"+name+"Oper"+num2str(quant)
string sCard="Clu_Car9_"+name+"Oper"+num2str(quant)
string sCenter="Clu_Cen9_"+name+"Oper"+num2str(quant)
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
Duplicate/FREE/O chargebag, cropy,wy
Duplicate/FREE/O chargebag, cropx,wx
make/FREE/O/D/N=(numpnts(chargebag),1) obs
Make/FREE/O/D/N=(Order) compcoef=1
cropy=Ciblebag[x][3]
cropx=ciblebag[x][1]
wy=Ciblebag[x][3]
wx=Ciblebag[x][1]
wave w_coef
variable k,n=100,test=1,t=ticks
for(k=0;k<n && test;k+=1)
	CurveFit/N/Q/M=2/W=0 poly Order, cropy/X=cropx
	compcoef=abs(w_coef-compcoef)
	test=wavemin(compcoef)
	compcoef=w_coef
	obs=abs(wy-poly(w_coef,wx))
	Cluster(nbClust,obs,3)
	wave voronass
	extract/O wy, cropy, voronass==0
	extract/O wx, cropx, voronass==0
endfor
Duplicate/O Means $sCenter
Duplicate/O Card $sCard
Duplicate/O Voronass $sClass
Killwaves Means,Card,Voronass
end

/////////////////////////////////////////////////
//K-Means Lloyd
/////////////////////////////////////////////////
function initiateKmeans(k,obs,init)
variable k,init
wave obs
variable dim=dimsize(obs,1)
variable n=dimsize(obs,0)
make/O/D/N=(k,dim) means
make/FREE/O/D/N=(n) trans
variable contdim
	switch(init)	// numeric switch
		case 1:
		for(contdim=0;contdim<dim;contdim+=1)
			trans=obs[x][contdim]
			means[][contdim]=mean(trans)+enoise(wavemax(trans)-wavemin(trans))
		endfor
			break
		case 2:
		make/O/D/N=(n) voronass=floor(abs(enoise(k)))
		newMeans(voronass,obs)
			break
		case 3:
		means=x
			break
	endswitch
end

function voronoi(means,obs)
wave means,obs
variable n=dimsize(obs,0)
variable k=dimsize(means,0)
variable dim=dimsize(obs,1)
make/FREE/O/D/N=(k) vorodist
make/O/D/N=(n) voronass
variable contn,contk,contdim
variable distance
for(contn=0;contn<n;contn+=1)
	for(contk=0;contk<k;contk+=1)
		distance=0
		for(contdim=0;contdim<dim;contdim+=1)
			distance+=(obs[contn][contdim]-means[contk][contdim])^2
		endfor
		vorodist[contk]=sqrt(distance)
	endfor
	findvalue/V=(wavemin(vorodist)) vorodist
	voronass[contn]=V_value
endfor
end

function newMeans(voronass,obs)
wave voronass,obs
variable n=dimsize(obs,0)
variable k=wavemax(voronass)+1
variable dim=dimsize(obs,1)
make/O/D/N=(k,dim) means=0
make/O/D/N=(k) card=0
variable contn,contdim
for(contn=0;contn<n;contn+=1)
	card[voronass[contn]]+=1
	for(contdim=0;contdim<dim;contdim+=1)
		means[voronass[contn]][contdim]+=obs[contn][contdim]
	endfor
endfor
means/=card[x]
end

function Cluster(k,obs,init)
variable k,init
wave obs
//newmovie/A/F=60
initiateKmeans(k,obs,init)
wave means, voronass
variable crit=1,iter
for(iter=0;crit!=0;iter+=1)
	duplicate/FREE/O means transmeans
	voronoi(means,obs)
	newMeans(voronass,obs)
	transmeans=transmeans^2-means^2
	crit=wavemax(transmeans)
	//doupdate
	//addmovieframe
endfor
make/O/D/N=(dimsize(means,0)) transvoronas=means[x][0]
wavetransform zapnans transvoronas
//print numpnts(transvoronas)
killwaves transvoronas
end
/////////////////////////////////////////////////

Function Calculate(mass,int,charj)
variable mass,int,charj
nvar charge, inco, nbprop, treedepth
charge=charj
inco=mass
wave/Z list_molperm, molecule
wave incox,incoy,incoxdm
incox[0]=inco
incoxdm[0]=defmass(inco)
incoy[0]=int
string lawave
analysemass("arg")
setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
lawave=stringfromlist(0,wavelist("molprop_*",";",""))
molecule=0
tuelesspec()
duplicate/O $lawave molecule
genestringmol()
ListBox list0 selRow=0, win=attributeur
controlinfo/W=dmvm popup0
generatetree(S_value,treedepth,inco)
end

function RemoveClassif(bagname,classif,criter)
string bagname
wave classif
variable criter
extract/FREE/INDX classif,mol2kill,classif==criter
extract/O classif,classif,classif!=criter
MolbagWithdrawFrom(mol2kill,bagname)
end

function MolbagDuplicate(name,newbagname)
string name,newbagname
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
string sciblenew="Ciblebag_"+newbagname
wave ciblebagnew=$sciblenew
variable n
if(waveexists(ciblebagnew))
	n=dimsize(ciblebagnew,0)
	make/FREE/O/D/N=(n) mol2kill=x
	MolbagWithdrawFrom(mol2kill,newbagname)
endif
Duplicate/FREE chargebag, index
index=x
MolbagExtractor(name,index,newbagname)
end

Window Somogyi_style_plot() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(339,203,709.5,674) roi2 vs roi0
	AppendToGraph/L=intensity roi1 vs roi0
	AppendToGraph/L=intensity treeintens vs tree
	AppendToGraph defmasstree vs tree
	AppendToGraph/L=intensity wave1stripped vs wave0stripped
	AppendToGraph wave2stripped vs wave0stripped
	ModifyGraph gbRGB=(34816,34816,34816)
	ModifyGraph mode(roi2)=3,mode(roi1)=8,mode(wave1stripped)=1,mode(wave2stripped)=3
	ModifyGraph marker(roi2)=19,marker(roi1)=19,marker(treeintens)=19,marker(defmasstree)=19
	ModifyGraph marker(wave2stripped)=19
	ModifyGraph lSize(treeintens)=2,lSize(defmasstree)=2,lSize(wave1stripped)=1.2
	ModifyGraph rgb(treeintens)=(65280,0,52224),rgb(defmasstree)=(65280,0,52224)
	ModifyGraph msize(roi2)=2,msize(roi1)=2
	ModifyGraph useMrkStrokeRGB(roi2)=1,useMrkStrokeRGB(roi1)=1,useMrkStrokeRGB(treeintens)=1
	ModifyGraph useMrkStrokeRGB(defmasstree)=1,useMrkStrokeRGB(wave2stripped)=1
	ModifyGraph zmrkSize(wave2stripped)={logwave1,*,*,1,5}
	ModifyGraph zColor(roi2)={roi1,-0.005,0.005,RedWhiteBlue256},zColor(roi1)={roi1,-0.005,0.005,RedWhiteBlue}
	ModifyGraph zColor(wave1stripped)={logwave1,*,*,YellowHot},zColor(wave2stripped)={logwave1,*,*,YellowHot}
	ModifyGraph log(intensity)=1
	ModifyGraph zero(intensity)=1
	ModifyGraph mirror=1
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	ModifyGraph lblPos(left)=69,lblPos(intensity)=69
	ModifyGraph lblLatPos(intensity)=5
	ModifyGraph freePos(intensity)=0
	ModifyGraph axisEnab(left)={0,0.48}
	ModifyGraph axisEnab(intensity)={0.52,1}
	Label left "Mass defect"
	Label bottom "Mass/Charge"
	Label intensity "Intensity"
	Cursor/P A roi2 -1
	ShowInfo
	Legend/C/N=text0/J/A=MC/X=18.51/Y=-17.98 "H\\S+\\M + (NH) C\\B\\f04\\f05X\\M\\f02 \\K(0,0,65280)(CH\\B2\\M)\\Bn\r\\s(wave1stripped) wave1stripped"
	AppendText "\\s(wave2stripped) wave2stripped"
	ColorScale/C/N=text1/F=0/B=1/A=MC/X=38.21/Y=33.39 trace=roi2, heightPct=40
	ColorScale/C/N=text1 tickLblRot=90, fstyle=1
	ColorScale/C/N=text1 userTicks={colorsappval,colorsapplab}
	ColorScale/C/N=text1 axisRange={-0.006,0.006,0}
EndMacro

function storeLayout()
string leswin=winlist("*",";","VISIBLE:1")
string lawin
variable k,n=itemsinlist(leswin)
Make/O/T/N=(n,5) WinLayout
for(k=0;k<n;k+=1)
	lawin=stringfromlist(k,leswin)
	getwindow $lawin, wsize
	WinLayout[k][]={{lawin},{num2str(V_left)},{num2str(V_top)},{num2str(V_right)},{num2str(V_bottom)}}
endfor
end

function setLayaout(winlayout)
wave/T winlayout
variable k,n=dimsize(winlayout,0)
string lawin
for(k=0;k<n;k+=1)
	lawin=winlayout[k][0]
	if(wintype(lawin))
		Movewindow/W=$lawin str2num(winlayout[k][1]),str2num(winlayout[k][2]),str2num(winlayout[k][3]),str2num(winlayout[k][4])
	endif
endfor
end

Function Mendev_Formula_Add(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
wave molecule
nvar charge
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			//Duplicate/O str2mol(sval) molecule
			molecule+=str2mol(sval)[x]
			charge = str2charge(sval)
			tuelesspec()
			genestringmol()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End


function/wave str2mol(str)
string str
variable k,n=strlen(str)
string elt1=" ", elt2
Make/FREE/O/N=(114) res=0
for(k=0;strlen(elt1) && k<(n+1);k+=1)
	Splitstring/E="([[:upper:]][[:lower:]]*)([[:punct:]]?[[:digit:]]*)" str, elt1,elt2
	str=str[strlen(elt1)+strlen(elt2),inf]
	if(strlen(elt1) && strlen(elt2))
		res[elt2indx(elt1)]+=str2num(elt2)
	elseif(strlen(elt1))
		res[elt2indx(elt1)]+=1
	endif
endfor
if(k>n)
	print "Bad format"
endif
return res
end

function str2charge(str)
string str
variable k,n=strlen(str),res=0
string elt1, elt2,elt3
for(k=0;strlen(str) && k<(n+1);k+=1)
	Splitstring/E="[[:alnum:]]*([[:space:]]+[[:punct:]]+[[:digit:]]*)" str, elt1
	if(strlen(elt1))
		if(stringmatch(elt1," +"))
			res+=1
		elseif(stringmatch(elt1," -"))
			res-=1
		else
			res+=str2num(elt1)
		endif
	endif
	str=str[strsearch(str,elt1,0)+strlen(elt1),inf]
endfor
return res
end

function elt2indx(elt)
string elt
wave elements
return whichlistitem(elt,wave2list(elements),";",0,0)
end

function/wave SplitDecimal(var)
variable var
string str=num2str(var),elt1,elt2
splitstring/E="([[:digit:]]+).([[:digit:]]+)" str,elt1,elt2
make/FREE/O/D/N=2 res
print elt1, elt2
res[0]=str2num(elt1)
res[1]=str2num(elt2)
return res
end

//////////////////////////////////////////////////////
///////// CHROMATOGRAMS
//////////////////////////////////////////////////////

////// LOADING CHROMATOGRAMS

function loadchroma()
string name="New_chromatogram"
prompt name, "Name this Sample"
doprompt "Enter a good name for this sample",name
variable t=ticks
if(V_flag==1)
	return -1
elseif(V_flag==0)
	string cheminversfichier
	string fichiery, fichierx, fichierz,fichiert
	GetFileFolderInfo/Q
	if(V_flag==-1)
		return -1
	elseif(V_flag==0)
		string nomfichier
		NewPath/O/Q ledossier ParseFilePath(1, S_Path, ":", 1, 0)
		nomfichier="chromafield1.itx"
		LoadWave/O/Q/T/D/P=ledossier nomfichier
		nomfichier="chromato1.txt"
		LoadWave/O/A=Chroma/Q/G/D/P=ledossier nomfichier
		fichiery=name[0,20]+"_1chroma"
		fichierx=name[0,20]+"_0chroma"
		fichierz=name[0,20]+"_2chroma"
		fichiert=name[0,20]+"_3chroma"
		Rename chromafield0, $fichierx
		Rename chromafield1, $fichiery
		Rename Chroma1, $fichierz
		Rename Chroma2, $fichiert
		wave Chroma0
		killwaves Chroma0
		display/N=TIC/K=1 $fichierz vs $fichiert
		if(waveexists(MaxTimeRange)==0)
			make/O/D/N=2 MaxTimeRange={wavemin($fichiert),wavemax($fichiert)}
		endif
		listelesdonnees()
	endif
endif
print "Loading time: " + num2str((ticks-t)/60)
end

////// BUILDING ION MAP

function ConstruitMap(name)
string name
string sChroma0=name+"_0chroma"
string sChroma1=name+"_1chroma"
string sChroma2=name+"_2chroma"
string sChroma3=name+"_3chroma"
string sChroma4=name+"_4chroma"
string sChroma5=name+"_5chroma"
string sChroma7=name+"_7chroma"
string sChroma13=name+"_13chroma"
Duplicate/FREE $sChroma0 chromafield0
Duplicate/FREE $sChroma1 chromafield1
Duplicate/FREE $sChroma2 tic
Duplicate/FREE $sChroma3 chromatime
variable t=ticks

// On dégomme la première colonne qui est une colonne idiote
deletepoints/M=1 0, 1, chromafield0, chromafield1
deletepoints 0,1, chromatime,tic

// On coupe en temps
wave MaxTimeRange //generated lors du chargement
duplicate/O MaxTimeRange indexTimeRangeMax
findlevel/P/Q Chromatime MaxTimeRange[0]
indexTimeRangeMax[0]=round(V_levelX)
findlevel/P/Q Chromatime MaxTimeRange[1]
indexTimeRangeMax[1]=round(V_levelX)
deletepoints/M=1 indexTimeRangeMax[1],dimsize(chromafield0,1), chromafield0
deletepoints/M=1 indexTimeRangeMax[1],dimsize(chromafield1,1), chromafield1
deletepoints/M=1 0,indexTimeRangeMax[0], chromafield0
deletepoints/M=1 0,indexTimeRangeMax[0], chromafield1
Extract/O Chromatime, Chromatime, x>=indexTimeRangeMax[0] && x<=indexTimeRangeMax[1]

variable msup=wavemax(chromafield0)
variable minf=wavemin(chromafield0)

// FAT Noise
// On détermine les niveaux de bruits
wave noizelevelz=fieldnoiselvl(chromafield1)
duplicate/O noizelevelz bruit666
extract/O noizelevelz, bruit666, noizelevelz>0
wave NoiseLevelAjust
variable maxnoise=10^(mean(bruit666)+NoiseLevelAjust[0]*sqrt(variance(bruit666)))
Print "Noise level : ",maxnoise," ("+num2str(NoiseLevelAjust[0])+" sigma)"
// On vérifie qu'on a un niveau de bruit qui existe
//if(numtype(maxnoise)!=0)
	//Abort("Error - Incorect noise level")
//else
	string newnoises
	variable tp=ticks
	prompt newnoises, "The current noise level is "+num2str(maxnoise)+" Click cancel to continue and use this value. Enter your custom noise level for this sample is you want to."
	Doprompt "Custom Noise Level",newnoises
	tp=ticks-tp
	if(V_Flag==0)
		variable newnoise=str2num(newnoises)
		if(numtype(newnoise!=0))
			Abort("Error - The new noise must be a number")
		else
			maxnoise=newnoise
			newnoises="Custom"
			Print "Custom noise level : ",maxnoise
		endif
	endif
//endif
// On multiplie les int & masses inférieures au seuil par 0. Sinon, mulitplie par 1
duplicate/FREE/O chromafield0 noizefield0
duplicate/FREE/O chromafield1 noizefield1 selector
Multithread selector=(chromafield1>maxnoise)
FastOP noizefield0=selector*noizefield0
FastOP noizefield1=selector*noizefield1
//if(sum(noizefield0)==0 || sum(noizefield1)==0 )
	//Abort("Error - Noise level too high: no data remaining")
//endif
// On normalise par rapport au seuil de bruit :
// Seuil de bruit = 1
//Multithread noizelevelz=10^noizelevelz
//Multithread chromafield1/=noizelevelz
//killwaves noizelevelz

// Détermination du numéro de spectre contenant le plus de points et traitement
variable numspectre=Look4Max(chromafield0, chromafield1, noizefield0, noizefield1)
diffmass4Max(chromafield0,numspectre)
// Rebinnage
wave diffmass, spectreMass
extractDM(chromafield0, diffmass, spectreMass,numspectre)
killwaves diffmass, spectreMass
// Construction de la map
wave Coef_Diffmass_mass,MaxTimeRange
ReconstructionMass(Coef_Diffmass_mass, chromafield0, noizefield0, noizefield1,minf,msup,chromatime)
wave combled
print "Map size (RT,m/z): ",dimsize(combled,1),dimsize(combled,0)
MatrixOP/O $sChroma7=sumrows(combled)
Duplicate/O Chromamass $sChroma4
Duplicate/O combled $sChroma5
Setscale y,wavemin(MaxTimeRange),wavemax(MaxTimeRange), $sChroma5
make/O/D/N=(6) $sChroma13
wave Chroma13=$sChroma13
Chroma13[0]=wavemin(MaxTimeRange) //min time
Chroma13[1]=wavemax(MaxTimeRange) //maxtime
Chroma13[2]=indexTimeRangeMax[0] // min time index
Chroma13[3]=indexTimeRangeMax[1] // max time index
Chroma13[4]=NaN // noise level of the map
Chroma13[4]=NoiseLevelAjust[0]*(stringmatch(newnoises,"Custom")!=1) // noise level of the map 
Chroma13[5]=maxnoise // noise cut
//killwaves combled, chromamass,Coef_Diffmass_mass
t=ticks-t
print "Building of the Ion map in "+num2str((t-tp)/60)+" secondes"
end

////// ION MAP FUNCTIONS

function/WAVE fieldnoiselvl(field1)
wave field1
variable k,n=dimsize(field1,1), m=dimsize(field1,0)
make/FREE/O/D/N=(n) noizelevelz
make/FREE/O/D/N=(m) transint
for(k=0;k<n;k+=1)
	transint=field1[x][k]
	noizelevelz[k]=seuilbruitFAT(transint)
endfor
return noizelevelz
end

function Look4Max(chromafield0, chromafield1, noizefield0, noizefield1)
wave chromafield0, chromafield1, noizefield0, noizefield1
// On boucle sur le chromafield en partant de la dernière ligne
// On boucle sur toutes les colonnes k de de la ligne l avant de remonter d'une ligne
// On s'arrête quand on a trouvé la première masse diffmass non nulle
variable col=dimsize(chromafield0,1), ligne=dimsize(chromafield0,0), k, l
variable dm, numspectre=inf
for(l=ligne-1;l<0;l-=1)
	for(k=col-1;k<0;k-=1)
		dm=abs(chromafield0[l-1][k]-chromafield0[l][k])
		// A-t-on trouvé le spectre avec le plus de points?
		if(dm!=0)
			// On enregistre la colonne du spectre avec le plus de points
			numspectre=k
		endif
	endfor
endfor
// On dégomme les lignes inutiles
Deletepoints/M=0 l,ligne, chromafield0, noizefield0, chromafield1, noizefield1
return numspectre
end

function diffmass4Max(chromafield0, numspectre)
wave chromafield0
variable numspectre
// On construit le DM pour le spectre numéro numspectre
// On extrait le spectre en question
variable taille=(dimsize(chromafield0,0))
make/O/D/N=(taille) spectreMass
Multithread spectreMass=chromafield0[x][numspectre]
// On calcule le DM
duplicate/O spectreMass diffmass
Multithread diffmass= abs(spectreMass[min(x+1,taille)]-spectreMass[x]) // erreur d'index ici (out of range)
// on enlève les zéros éventuels
duplicate/FREE/O diffmass transdiffmass
extract/O diffmass, diffmass, transdiffmass>0
extract/O spectreMass, spectreMass, transdiffmass>0
end

function extractDM(chromafield0, diffmass, spectreMass,numspectre)
wave chromafield0, diffmass, spectreMass
variable numspectre
// On récupère min et max du spectre
make/FREE/O/D/N=(2) DMMinMax
DMMinMax[0]=wavemin(spectreMass)
DMMinMax[1]=wavemax(spectreMass)
variable p=10 // Pas de masse pour crop le DM
variable k, l, n=round((DMMinMax[1]-DMMinMax[0])/p-1)
variable mass=DMMinMax[0], massindex=0
make/O/D/FREE/N=(n) cropmass=Nan, cropdiffmass=Nan
for(l=0;l<numpnts(spectreMass);l+=1) // On screen toutes les masses
	// La différence entre les deux masse est supérieure au pas choisi?
	if(abs(mass-spectreMass(l))>=p)		
		// On sélectionne la plus petite valeur de diffmass entre les deux index
		cropdiffmass[k]=wavemin(diffmass, massindex,l)
		findlevel/P/Q/R=(massindex,l) diffmass cropdiffmass[k]
		cropmass[k]=spectreMass[V_levelX]
		// On réassigne les masses
		mass=spectreMass(l)
		massindex=l
		// On incrémente le compteur du crop
		k+=1
	endif
endfor
// On arrivera surement pas à remplir "tous les pas" : il reste des vides
wavetransform zapNaNs cropmass
wavetransform zapNans cropdiffmass
// On fit l'enveloppe trouvée
make/O/D/N=3 w_coef=0
CurveFit/H="110"/Q/M=2/N/NTHR=0/W=0 poly 3, kwCWave=w_coef ,cropdiffmass/X=cropmass
duplicate/O w_coef Coef_Diffmass_mass
end

function ReconstructionMass(Coef_Diffmass_mass, chromafield0, noizefield0, noizefield1,minf,msup,chromatime)
wave Coef_Diffmass_mass, chromafield0, noizefield0, noizefield1,chromatime
variable minf,msup

// Construction de chromafield2
duplicate/FREE/O chromafield0 chromafield2
Multithread chromafield0=chromafield0/chromafield0*chromafield0
minf=wavemin(chromafield0)
msup=wavemax(chromafield0)
Multithread chromafield2=round(-1/(Coef_Diffmass_mass[2]*noizefield0)+1/(Coef_Diffmass_mass[2]*minf))

wave indexTimeRangeMax
// Construction du field des intensités dans le champ des indices
make/FREE/O/D/N=(numpnts(chromafield2)) chromaline
Multithread chromaline=chromafield2
// On dégomme les valeurs négatives
extract/O chromaline, chromaline, chromaline>0
// On les range dans l'ordre
sort chromaline, chromaline
// On dégomme les valeurs identiques
duplicate/FREE/O chromaline transchromaline
Multithread transchromaline=abs(chromaline[x]-chromaline[x+1])
extract/O chromaline, chromaline, transchromaline>0 
// On ajoute des lignes là où c'est nécessaire
duplicate/O/FREE chromaline transchromaline, transpos
Multithread transchromaline=abs(chromaline[x]-chromaline[x+1])
extract/O/INDX chromaline, transpos, transchromaline>1 
variable posl
variable Vvalue, Vindex
for(posl=numpnts(transpos)-1;posl>0;posl-=1)
	Vvalue=chromaline[transpos[posl]]+1
	Vindex=transpos[posl]
	insertpoints Vindex+1,1,chromaline
	chromaline[Vindex+1]=Vvalue
endfor
// On converti les index uniques en leurs masses
duplicate/O chromaline chromamass
multithread chromamass=1/(1/minf-Coef_Diffmass_mass[2]*chromaline)

// Maintenant qu'on a les indices uniques, on build les intensités en field!
// On check chaque indice dans chromafield2, on détermine sa position (ligne) dans recint grâce à chromaline et on somme son intensité
// On précalcule les positions des index dans recint
Make/FREE/O/D/N=(wavemax(chromaline)) pos=Nan
variable k, l=numpnts(chromaline)
for(k=0;k<l;k+=1)
	pos[chromaline[k]]=k
endfor
// On génère recint
if(numtype(dimsize(chromaline,0))!=0)
	Abort("Error - Trying to build a map with no or too many points")
endif
make/O/FREE/D/N=(dimsize(chromaline,0),abs(indexTimeRangeMax[1]-indexTimeRangeMax[0])) recint
multithread recint=0
recintgenvec(recint,pos,chromafield2,noizefield1)
// On comble les trous dus à l'échantillonage
comble(recint,chromaline)
end

threadsafe function recintgencalc(xx,yy,recint,pos,chromafield2,noizefield1)
variable xx,yy
wave recint,pos,chromafield2,noizefield1
recint[pos[chromafield2[xx][yy]]][yy]+=noizefield1[xx][yy]
return 0
end

function recintgenvec(recint,pos,chromafield2,noizefield1)
wave recint,pos,chromafield2,noizefield1
variable m=dimsize(recint,1),n=dimsize(chromafield2,0)
MAKE/O/D/N=(n,m)/FREE ghost
multithread ghost=recintgencalc(x,y,recint,pos,chromafield2,noizefield1)
end 

function comble(unemap,indexlist)
wave unemap, indexlist
duplicate/O unemap combled
// Max d'un point nul situé entre deux points en masse non nuls à un temps donné pour des index en masses continus
Multithread combled+=(indexlist[x+1]-indexlist[x]==1)*(max(unemap[x-1][y],unemap[x+1][y]))*((unemap[x+1][y]!=0 && unemap[x-1][y]!=0) && unemap[x][y]==0)
//Smooth/DIM=0/E=1/EVEN/B=1 3, combled
// Smooth en temps: Combinaison d'un smooth gaussien sur 5 points et d'une moyenne glissante sur 3 points
//Smooth/DIM=1/E=1 4, combled
//Smooth/DIM=1/E=1/EVEN/B 5, combled
//Smooth/E=1 3, combled
end

////// END ION MAP FUNCTIONS

////// ISLAND DETECTION - HOSHEN-KOPELMAN FUNCTIONS

Function ChromatoMap2vectorsHK(name)
string name
string sMap=name+"_5chroma"
string sHKRes=name+"_6chroma"
//string sHKParents=name+"_10chroma"
wave Map=$sMap
if(waveexists($sMap))
	Duplicate/O hoshenkopelman(name,Map) $sHKRes
	wave HKRes=$sHKRes
	insertpoints/M=1 10,1,HKRes
	HKRes[][9]=defmass(HKRes[x][6])
endif
end

////// PEAK DETECTION - HOSHEN-KOPELMAN FUNCTIONS

threadsafe function/WAVE hoshenkopelman(name,w)
string name
wave w
variable t=ticks
Duplicate/O/FREE w, hk
variable n=dimsize(w,0)
variable m=dimsize(w,1)
variable k,l,j,count
variable left, above
extract/O/INDX/FREE w, transx, w>0
duplicate/O/FREE transx, transy, labels,intensitouille
multithread transy=floor(transx/n)
multithread transx=transx-transy*n
multithread intensitouille=w[transx][transy]
insertpoints numpnts(labels),1,labels
count=dimsize(labels,0)
multithread labels=x
FastOP hk=(count)
for(k=0;k<count-1;k+=1)
	left=w[max(0,transx[k])][max(0,transy[k]-1)]//w[max(1,transx[k])][max(1,transy[k]-1)]
	above=w[max(0,transx[k]-1)][max(0,transy[k])]//w[max(1,transx[k]-1)][max(1,transy[k])]
	if(left==0 && above==0)
		hk[transx[k]][transy[k]]=k
	elseif(left!=0 && above==0)
		hk[transx[k]][transy[k]]=findHK(labels,hk[transx[k]][transy[k]-1])
		if(transy[k]==0)
			//print transx[k]
			hk[transx[k]][transy[k]]=k
		endif
	elseif(left==0 && above!=0)
		hk[transx[k]][transy[k]]=findHK(labels,hk[transx[k]-1][transy[k]])
	elseif(left!=0 && above!=0)
		if(transy[k]==0)
		else
			unionHK(labels,hk[transx[k]-1][transy[k]],hk[transx[k]][transy[k]-1])
		endif
		hk[transx[k]][transy[k]]=findHK(labels,hk[transx[k]-1][transy[k]])
	endif
endfor
Multithread hk=findHK(labels,hk)
make/FREE/O/D/N=(count,9) res
multithread res=0
make/FREE/O/D/N=(count,1) transint
multithread res[][2]=inf
multithread res[][3]=-inf
multithread res[][4]=inf
multithread res[][5]=-inf
multithread res[][7]=0
for(k=0;k<count;k+=1)
	res[hk[transx[k]][transy[k]]][0]+=1//pixel number
	res[hk[transx[k]][transy[k]]][1]+=w[transx[k]][transy[k]]//integrated intensity
	res[hk[transx[k]][transy[k]]][2]=min(res[(hk[transx[k]][transy[k]])][2],transx[k])//mass index
	res[hk[transx[k]][transy[k]]][3]=max(res[(hk[transx[k]][transy[k]])][3],transx[k])//mass index
	res[hk[transx[k]][transy[k]]][6]+=transx[k]//mass index
	res[hk[transx[k]][transy[k]]][4]=min(res[(hk[transx[k]][transy[k]])][4],transy[k])//time index
	res[hk[transx[k]][transy[k]]][5]=max(res[(hk[transx[k]][transy[k]])][5],transy[k])//time index
	if(w[transx[k]][transy[k]]>transint[hk[transx[k]][transy[k]]])
		res[hk[transx[k]][transy[k]]][7]=transy[k] //time index at max intensity
		transint[hk[transx[k]][transy[k]]]=w[transx[k]][transy[k]]// max intensity of the peak (not stored)
	endif
	res[hk[transx[k]][transy[k]]][8]=hk[transx[k]][transy[k]]//HKParent
endfor
Multithread res[][6]=res[x][6]/res[x][0]
ConvertIndx2Mass(name,res)
ConvertIndx2Time(name,res)
//*****Filtering*****//
variable setfiltering=1
// Deleting null values
extract/O/FREE/INDX res,izfull,res>0 && q==0
make/O/FREE/D/N=(dimsize(izfull,0),dimsize(res,1)) res2
multithread res2=res[izfull[x]][y]
if(setfiltering==1)
// Background parent
extract/O/FREE res2,transfilter,q==0
variable BackgroundParent=wavemax(transfilter)
findvalue/V=(BackgroundParent) transfilter
BackgroundParent=hk[V_value-(floor(V_value/dimsize(transfilter,0)))*dimsize(transfilter,0)][8]
// START FILTERS
// Filtering HK results based on some arbitrary parameters
// 1) If delta(m)=0u, delete peaks from results
Duplicate/O/FREE HKFilter(res2,hk,labels,BackgroundParent, 2,3,0, "="), res2
// 2) If delta(t)<10 time points, delete peaks from results
variable deltat=dimdelta(w,1)
Duplicate/O/FREE HKFilter(res2,hk,labels,BackgroundParent, 4,5,deltat*10, "<"), res2
// 3) If pixel number <30, delete peaks from results
// Parents des pics à retirer
Extract/O/FREE res2,transfilter,(res2[x][0])<30 && q==8
j=dimsize(transfilter,0)
for(l=0;l<j;l+=1)
	unionHK(labels,transfilter[l],BackgroundParent)
endfor
// Index des pics à conserver
Extract/O/FREE/INDX res2,transfilter,(res2[x][0])>=30 && q==0
Make/O/D/FREE/N=(dimsize(transfilter,0),dimsize(res2,1)) res
Multithread res=res2[transfilter[x]][y]
Duplicate/O/FREE res, res2
// END FILTERS
endif
Make/O/I/FREE/N=(count,4) leanmap
Multithread leanmap[][0]=transx[x]
Multithread leanmap[][1]=transy[x]
Multithread leanmap[][2]=hk[transx[x]][transy[x]]
Multithread leanmap[][3]=intensitouille[x]
string sleanmap=name+"_12chroma"
duplicate/O/I leanmap, $sleanmap
print "Island detection in ", (ticks-t)/60, " sec"
return res2
end

// *********
// Regénération multithread de la map
function mapgenerator(name,selector)
string name, selector
// Selector = "mass" or "island"
variable l,m,t=ticks
string stimewave=name+"_3chroma"
string smasswave=name+"_4chroma"
string ionmap=name+"_5chroma"
string hkmap=name+"_10chroma"
string sleanmap=name+"_12chroma"
string smapparam=name+"_13chroma"
wave timewave=$stimewave
wave masswave=$smasswave
wave leanmap=$sleanmap
wave mapparam=$smapparam
// Récupération de la dimension temporelle (index)
variable mintime=mapparam[2],maxtime=mapparam[3]
// Génération de la map
variable massdim=dimsize(masswave,0)
variable timedim=maxtime-mintime
variable leandim=dimsize(leanmap,0)
make/O/I/FREE/N=(massdim,timedim) mapgen
Multithread mapgen=0
if(stringmatch(selector,"mass"))
	vecmapgen(leanmap,mapgen,3)
	duplicate/O/I mapgen, $ionmap
	Setscale/I y,mapparam[0],mapparam[1],$ionmap
elseif(stringmatch(selector,"island"))
	vecmapgen(leanmap,mapgen,2)
	duplicate/O/I mapgen, $hkmap
	Setscale/I y,mapparam[0],mapparam[1],$hkmap
else
	print "Error - the selector must be 'mass' or 'island'"
endif
//print "Temps de regénération : ",(ticks-t)/60," sec"
end

threadsafe function leanvect(xx,adresses,dest,col)
variable xx,col
wave adresses,dest
dest[adresses[xx][0]][adresses[xx][1]]=adresses[xx][col]
return 0
end

function vecmapgen(adresses,dest,col)
wave adresses,dest
variable col
variable n=dimsize(adresses,0)
MAKE/O/D/N=(n)/FREE ghost
multithread ghost=leanvect(x,adresses,dest,col)
end
// END - Regénération multithread de la map
// *********

Threadsafe Function/Wave HKFilter(res2,hk,labels,BackgroundParent, col1,col2,threshold, mode)
wave res2, hk,labels
variable BackgroundParent, col1,col2,threshold
string mode
variable j,l
// Mode="=" : use the filter to check for values equal to the threshold value
// Mode=">" : use the filter to check for values higher than the threshold value
// Mode="<" : use the filter to check for values lower than the threshold value
// Parents des pics à retirer
if(cmpstr(mode,"=")==0 || cmpstr(mode,"<")==0  || cmpstr(mode,">")==0 )
	if(cmpstr(mode,"=")==0)
		extract/O/FREE res2,transfilter,(res2[x][col2]-res2[x][col1])==threshold && q==8
	elseif( cmpstr(mode,">")==0)
		extract/O/FREE res2,transfilter,(res2[x][col2]-res2[x][col1])>threshold && q==8
	elseif(cmpstr(mode,"<")==0)
		extract/O/FREE res2,transfilter,(res2[x][col2]-res2[x][col1])<threshold && q==8
	endif
	j=dimsize(transfilter,0)
	for(l=0;l<j;l+=1)
		unionHK(labels,transfilter[l],BackgroundParent)
	endfor
	// Index des pics à conserver
	if(cmpstr(mode,"=")==0)
		extract/O/FREE/INDX res2,transfilter,(res2[x][col2]-res2[x][col1])!=threshold && q==0
	elseif( cmpstr(mode,">")==0)
		extract/O/FREE/INDX res2,transfilter,(res2[x][col2]-res2[x][col1])<threshold && q==0
	elseif(cmpstr(mode,"<")==0)
		extract/O/FREE/INDX res2,transfilter,(res2[x][col2]-res2[x][col1])>threshold && q==0
	endif
	make/O/D/FREE/N=(dimsize(transfilter,0),dimsize(res2,1)) res
	multithread res=res2[transfilter[x]][y]
	duplicate/O/FREE res, res2

	return res2
else
	print "Error: not a valid mode"
	return res2 //no change on the initial matrix
endif
end

Threadsafe Function ConvertIndx2Mass(name,res)
string name
wave res
string sMasse=name+"_4chroma" 
wave mass=$sMasse
// Index2mass
multithread res[][2]=mass[res[x][2]]
multithread res[][3]=mass[res[x][3]]
multithread res[][6]=mass[res[x][6]]
end

Threadsafe Function ConvertIndx2Time(name,res)
string name
wave res
string sTime=name+"_3chroma" 
string sMap=name+"_5chroma" 
wave chromatime=$sTime
wave chromamap=$sMap 
// Time offset
variable timeoffset=dimoffset(chromamap,1)
				  
findlevel/Q chromatime, timeoffset
variable timeoffsetindx=ceil(V_levelX)
					  
// Index2Time
multithread res[][4]=chromatime[res[x][4]+timeoffsetindx]
multithread res[][5]=chromatime[res[x][5]+timeoffsetindx]
multithread res[][7]=chromatime[res[x][7]+timeoffsetindx]
end

threadsafe Function findHKp(parent,k)
Variable k
Wave parent	
	If(parent[k]!=k)
		parent[k]=findHKp(parent,parent[k])
	endif
	return parent[k]
end

threadsafe Function findHK(parent,k)
Variable k
Wave parent	
	for(k=k;parent[k]!=k;k=k)
		k=parent[k]
	endfor
	return parent[k]
end

threadsafe function unionHK(parent,pos1,pos2)
wave parent
variable pos1, pos2
parent[findHK(parent,pos1)]=findHK(parent,pos2)
end

// Classification KMeans sur les résultats de Hoshen-Kopelman
function classilots(hkres,n,displaygraph)
variable n, displaygraph
wave hkres
variable k=dimsize(hkres,0)
make/O/D/N=(k,5) crophk
crophk[][0]=hkres[x][0]
//crophk[][1]=hkres[x][1]
crophk[][1]=hkres[x][3]-hkres[x][2]
crophk[][2]=hkres[x][5]-hkres[x][4]
crophk[][3]=hkres[x][9]
matrixOP/FREE pop= crophk^t
variable dim=dimsize(pop,0)
make/O/D/N=(dim,n)/FREE classes
Kmeans/OUT=2/INW=classes pop
//edit/K=1 M_KMClasses
//edit/K=1 W_KMMembers
if(displaygraph==1)
display/K=1 hkres[][7] vs hkres[][6]
ModifyGraph mode=3,marker=19,useMrkStrokeRGB=1;DelayUpdate
display/K=1 hkres[][1] vs hkres[][6]
ModifyGraph log(left)=1
ModifyGraph mode=3,marker=19,useMrkStrokeRGB=1;DelayUpdate
endif
end

// Extraction des ilots détectés par HK selon les formules stoechiométriques renseignées dans un molbag
function extractPeakfromHK(shkname, sbagname)
string shkname, sbagname
wave hk=$shkname+"_6chroma"
wave Ciblebag=$"Ciblebag_"+sbagname
string sChromabag="Chromabag_"+sbagname
variable k=0,l,n=dimsize(Ciblebag,0),m=dimsize(hk,0),p
Duplicate/O/D hk, $sChromabag
wave Chromabag=$sChromabag
// Calcul des masses exactes théoriques
make/O/N=(n,1) transmass
multithread transmass=Ciblebag[x][1]/(1-(Ciblebag[x][3]/10^6))
// On recherche les pics chromato dont leur étendue en masse comprend la masse théorique recherchée
for(k=m;k!=-1;k-=1)
	l=0
	p=0
	for(l=0;l<n;l+=1)
		if(transmass[l]<hk[k][3] && transmass[l]>hk[k][2])
			// Pic à conserver
			p=1
			Chromabag[k][9]=transmass[l]
			// On sort de la boucle
			l=n
		endif
	endfor
	if(p==0)
		deletepoints/M=0 k,1,Chromabag
	endif
endfor
end

////// END ISLAND DETECTION - HOSHEN-KOPELMAN FUNCTIONS

////// PEAK DETECTION - PERSISTENT HOMOLOGY FUNCTIONS

// Découpe une portion d'ionmap autour des min/max masse et temps pour un ilotID donné
function island2chroma(ilotID, samplename)
variable ilotID
string samplename
wave time2indx=$samplename + "_3chroma"//temps <-> time index
wave mass2indx=$samplename + "_4chroma" //mass <-> mass index
wave ionmap=$samplename + "_5chroma" //ion map
wave hkresults=$samplename + "_6chroma"//HK table
wave hkmap=$samplename + "_10chroma" //HK map
wave leanmap=$samplename + "_12chroma"//lean map
wave mapgenparam=$samplename + "_13chroma"//map generation parameters


//variable timedim=dimsize(ionmap,1)
//variable massdim=dimsize(ionmap,0)
variable islandnum=dimsize(hkresults,0)
make/O/D/FREE/N=(islandnum) transnum
multithread transnum=hkresults[x][8]
// on trouve la ligne dans HK results de l'ilot en question
findlevel/Q transnum, ilotID //**//
if(V_flag==1)
	Abort("Error - Island not found\r\rIsland not found error could be due to:\r- Dirac in mass;\r- Less than 10 points in time;\r- Too large mass range;\r- The background.")
else 
	variable ilotline=round(V_levelX)
endif

wave meanmass
// Offset de la map en temps
variable mapoffset=mapgenparam[0] //min time of the map
variable minmass = hkresults[ilotline][2]
variable maxmass = hkresults[ilotline][3]
variable mintime = hkresults[ilotline][4]-mapoffset
variable maxtime = hkresults[ilotline][5]-mapoffset
meanmass[0]=hkresults[ilotline][6]
//print minmass, maxmass,mintime,maxtime

minmass=findmassortimevalue(minmass,mass2indx)//*min mass indx*//
maxmass=findmassortimevalue(maxmass,mass2indx)//*max mass indx*//
mintime=findmassortimevalue(mintime,time2indx)//*min time indx*//
maxtime=findmassortimevalue(maxtime,time2indx)//*max time indx*//

// Islandchroma = vecteur 1D qui somme les intensités pour différentes masses appartenant à l'ilotID
make/O/D/N=(maxtime-mintime) islandchroma=0
variable l,m
for(l=minmass;l<=maxmass;l+=1)
	for(m=mintime;m<=maxtime;m+=1)
		islandchroma[m-mintime]+=ionmap[l][m]*(hkmap[l][m]==ilotID)
	endfor
endfor

setscale x,hkresults[ilotline][4],hkresults[ilotline][5],islandchroma
end

function findmassortimevalue(val,w)
variable val
wave w
findlevel/Q w, val
if(V_flag==1)
	val=0
else 
	val=round(V_levelX)
endif
return val
end


// Algorithme Persistent Homology
// More informations :
// https://www.sthu.org/blog/13-perstopology-peakdetection/index.html
function pershomo(wi)
wave wi
duplicate/O/FREE wi, w, temps
temps=x
// On décime le nombre de points par un lissage agressif
// Pour les ilots peu étendus, paramètre minimum de 3
smooth (max(floor(numpnts(wi)/20),3)), w
duplicate/O w pli
// pli = 1, 0 ou -1
// 1 => Maximum local
// 0 => Valeur transitoire
// -1 => Minimum local
pli=((w[p-1]<=w[p])&&(w[p+1]<=w[p]))-((w[p-1]>=w[p])&&(w[p+1]>=w[p]))
extract/O/INDX pli, extrindx, pli!=0
extract/O temps, extrtime, pli!=0
extract/O pli, extrpli, pli!=0
extract/O w, extr, pli!=0
variable n=numpnts(extr),k
Duplicate/O extr,sortRindx
makeindex/R extr, sortRindx
Duplicate/O extr,birth,death,persi
// Birth et death prennent les valeurs d'intensités du point
// On attribue à chaque maximum local sa valeur d'intensité (Birth)
// Un minimum local est forcément encadré de deux maximas locaux
// La valeur du minimim local est attribuée au plus petit des minimums locaux (death)
birth=0
death=0
variable lepoint
for(k=0;k<n;k+=1)
lepoint=sortRindx[k]
	if(extrpli[lepoint]==1)
		birth[lepoint]=extr[lepoint]
	else
		if(birth[lepoint+1]>birth[lepoint-1])
			death[lepoint-1]=extr[lepoint]
		else
			death[lepoint+1]=extr[lepoint]
		endif
	endif
endfor
// Une persistance élevée (persi) indique un point qui a peu subit de death, donc un sommet de pic
persi=birth-death
// Correction des effets de bords
variable dimpersi=numpnts(persi)
persi[0]=0
persi[dimpersi]=0
if(persi[dimpersi-1]==0 && persi[dimpersi-2]==0)
	deletepoints/M=0 dimpersi-1,1,persi,extrtime,death,birth,extr
endif
if(persi[0]==0 && persi[1]==0)
	deletepoints/M=0 0,1,persi,extrtime,death,birth,extr
endif
// Nettoyage des persistances non nulles très faibles quand on a au moins 2 pics détectés
variable sizepersi=numpnts(persi)
//print "Detected peaks before filtering: ",(numpnts(persi)-1)/2
if(sizepersi>1)
	k=sizepersi-2
	for(k=sizepersi-2;k>0;k-=2)
		if(birth[k]<2*death[k])
			deletepoints/M=0 k,2,persi,extrtime,death,birth,extr
		endif
	endfor
endif

//print "Detected peaks after filtering: ",(numpnts(persi)-1)/2

// représentations : ilsandchroma vs extrtime (line) + persi vs extrtime (stick to zero + mark)
end

////// END PEAK DETECTION - PERSISTENT HOMOLOGY

////// PEAK FITING FUNCTIONS

// Fonction qui découpe les ilots, en fait des chromato, détecte les pics, les modélise et créé une wave de résultats
function island2chromastruct(samplename)
string samplename
wave hkresults=$samplename + "_6chroma"
variable ilotcount=dimsize(hkresults,0)
variable mintime,maxtime
// Note : le fond n'est pas dans HKresults
// Génération de _10chroma
string hkmap=samplename+"_10chroma"
mapgenerator(samplename,"island")
variable k=0,l,t=ticks
variable ilotID,coefcount,modelcount=0
variable fitparam=4 //H,c,l,exp
make/O/D/N=(ilotcount*10,11) $samplename + "_11chroma"
wave fitresults=$samplename + "_11chroma"
make/O/D/N=1 meanmass
for(k=0;k<ilotcount;k+=1)
	ilotID=hkresults[k][8]
	// Découpe
	island2chroma(ilotID,samplename)
	if(wavemax(islandchroma)/wavemin(islandchroma)>3)
		// Détection
		pershomo(islandchroma)
		// Génération des guess
		Fitguessfrompershomo(persi,extrtime,0)
		// On s'assure qu'on a bien détecté un pic
		if(numpnts(persi)>1)
			// Fit
			fitisland(islandchroma, fitinit)
			// Résultats
			l=0
			wave fitinit,fitcoef
			mintime=pnt2x(islandchroma,0)
			maxtime=pnt2x(islandchroma,numpnts(islandchroma))
			coefcount=numpnts(fitcoef)
			for(l=0;l<coefcount;l+=fitparam)
				fitresults[modelcount][0]=ilotID //IlotID
				fitresults[modelcount][1]=meanmass[0] //mean mass of the island
				// Paramètres après fit
				fitresults[modelcount][2]= fitcoef[l]*(fitcoef[l]>0)//H
				fitresults[modelcount][3]= abs(fitcoef[l+1])*(fitcoef[l+1]>mintime && fitcoef[l+1]<maxtime)//c
				fitresults[modelcount][4]= abs(fitcoef[l+2])//l
				fitresults[modelcount][5]= fitcoef[l+3]//exp
				// Valeurs d'initialisation
				fitresults[modelcount][6]= fitinit[l]//H
				fitresults[modelcount][7]= fitinit[l+1]//c
				fitresults[modelcount][8]= abs(fitinit[l+2])//l
				fitresults[modelcount][9]= fitinit[l+3]//exp
				// Calcul des aires des pics fittés
				fitresults[modelcount][10]=fitarea4island(samplename,islandchroma,modelcount)
				// Passage à la ligne suivante
				modelcount+=1
			endfor
		endif
	endif
endfor
deletepoints/M=0 modelcount, ilotcount*10-modelcount, fitresults
// Génération des positions temporelles des fits pour l'ion map
islandfit4ionmap(samplename)
print "Peaks fiting in ", (ticks-t)/60," sec"
killwaves $hkmap
end

// Fonction qui permet de fitter plusieurs gaussienes
// Prend une fonction 1D en entrée pour les coeff
function polygau(w,x) : fitfunc
wave w//vecteur obligatoire H,l,c
variable x
variable n=dimsize(w,0), k
variable res=0
//k+=3 => dépendant du nombre de paramètres
for(k=0;k<n;k+=3)
	// Somme de Gaussiennes
	//vecteur obligatoire H,l,c
	res+=w[k]*exp(-(x-w[k+2])^2/(2*w[k+1]^2))
endfor
return res
end

// Fonction qui permet de fitter plusieurs EMG
// Prend une fonction 1D en entrée pour les coeff
function polyEMG(w,x) : fitfunc
wave w//vecteur obligatoire H,c,l,exp
variable x
variable n=dimsize(w,0), k
variable res=0
//k+=4 => dépendant du nombre de paramètres
for(k=0;k<n;k+=4)
	// Somme de Exponentially modified Gaussian distribution
	// vecteur obligatoire H,c,l,exp
	res+=w[k]*exp(w[k+3]/2*(2*w[k+1]+w[k+3]*w[k+2]^2-2*x))*erfc((w[k+1]+w[k+3]*w[k+2]^2-x)/(sqrt(2)*w[k+2]))
endfor
return res
end

// Fonction qui permet de fitter plusieurs Pearson IV
// Prend une fonction 1D en entrée pour les coeff
function polyP4(w,x) : fitfunc
wave w//vecteur obligatoire H,c,l,forme,asym
variable x
variable n=dimsize(w,0), k
variable res=0
for(k=0;k<n;k+=5)
	// Somme de Pearson IV
	// vecteur obligatoire H,c,l,forme,asym
	res+= w[k]*(1+((x-w[k+1])/w[k+2])^2)^(-w[k+3])*exp(-w[k+4]*atan(((x-w[k+1])/w[k+2])))
endfor
return res
end

// Génération des guess à partir des données persistent homology
function Fitguessfrompershomo(persi,extrtime,mode)
wave persi, extrtime
variable mode //0=nettoyage, 1=pas de nettoyage
wave extr
// Nettoyage par amplitudes
if(mode==0)
	variable ampkiller=2 // PARAMETRE ARBITRAIRE
	extract/O/D/FREE extr,transextr,persi==0
	variable signalkiller=mean(transextr),sizepersi=numpnts(persi),k
	if(sizepersi>1)
		k=sizepersi-2
		for(k=sizepersi-2;k>0;k-=2)
			if(extr[k]<ampkiller*signalkiller)
				deletepoints/M=0 k,2,persi,extrtime,death,birth,extr
			endif
		endfor
	endif
endif
// On vérifie que l'on peut bien construire ce dont on a besoin
if(numpnts(persi)<=1)
	//Abort("Error - No peaks detected")
Else
	// Construction wave 1D des coefficients basé sur persi et extrtime
	variable n=numpnts(persi)
	variable fitparam=4 //H,c,l,exp
	// Nombre de coeff = nombre de paramètres de fit * nombre de fit
	variable coefnum=((n-1)/2)*fitparam
		make/O/D/N=(coefnum) fitinit
	variable fitcount=0
	for(k=0;k<coefnum;k+=fitparam)
		//H
		fitinit[k]=persi[2*fitcount+1]/2
		//c
		fitinit[k+1]=extrtime[2*fitcount+1]
		//l
		fitinit[k+2]=(extrtime[2*fitcount+2]-extrtime[2*fitcount])/10
		//exp
		fitinit[k+3]=1
		// Compteur de fit
		fitcount+=1
	endfor
endif
end

// Génération du fit
function fitisland(chroma,fitinit)
wave chroma, fitinit
duplicate/O fitinit, fitcoef
variable V_FitMaxIters=100
Funcfit/C/Q/N=1/G/NTHR=0/W=2 polyEMG fitcoef islandchroma/D
end

// Génération du guess initial
function fitislandguess(chroma,fitinit)
wave chroma, fitinit
duplicate/O fitinit, fitcoef
variable V_FitMaxIters=100
Funcfit/O/C/Q/N=1/G/NTHR=0/W=2 polyEMG fitcoef islandchroma/D
end

// Calcul de l'aire du pic chromato
function fitarea4island(name,chromato,fitline)
string name
wave chromato
variable fitline
string sfitresults=name + "_11chroma"
wave fitresults=$sfitresults
// Calcul de l'aire à partir des fits individuels
// Récupération des paramètres
make/O/D/FREE/N=(4) transcoef
transcoef[0]=fitresults[fitline][2]
transcoef[1]=fitresults[fitline][3]
transcoef[2]=fitresults[fitline][4]
transcoef[3]=fitresults[fitline][5]
// Création d'un chromatogramme simulé
duplicate/O/FREE islandchroma, transislandchroma
transislandchroma=0
transislandchroma=polyEMG(transcoef,x)
make/O/D/FREE/N=(dimsize(transislandchroma,0)) transtime
transtime=dimoffset(transislandchroma,0)+x*dimdelta(transislandchroma,0)
variable res=areaXY(transtime,transislandchroma)
return res
end

////// END PEAK FITING FUNCTIONS

////// 3D CHROMATOGRAMS

function ChromaD(name,mzcenter,mzwidth,stats)
string name
variable mzcenter, mzwidth, stats
string sTime=name+"_3chroma" 
string sMasse=name+"_4chroma"
string sMap=name+"_5chroma" 
duplicate/O/FREE $sTime chromatime 
duplicate/O/FREE $sMasse mass 
duplicate/O/FREE $sMap chromamap 
// Détection des bornes en index
variable minmass, maxmass
minmass=mzcenter-mzwidth/2
maxmass=mzcenter+mzwidth/2
findlevel/Q mass, minmass //**//
if(V_flag==1)
	variable minmindx=0
else 
	minmindx=round(V_levelX)
endif
findlevel/Q mass, maxmass //**//
if(V_flag==1)
	variable maxmindx=numpnts(mass)
else
	maxmindx=round(V_levelX)
endif
// On coupe la map en masse
deletepoints/M=0 maxmindx,dimsize(chromamap,0), chromamap, mass
deletepoints/M=0 0,minmindx, chromamap, mass
// Linéarisation des index de la map 
extract/O/FREE/INDX chromamap, trans, chromamap!=0
// Linéarisation de l'intensité 
extract/O/FREE chromamap, transz, chromamap!=0 
duplicate/O/FREE trans, transx, transy 
// Extraction des index de ligne et de colonne 
Multithread transy=floor(trans/dimsize(chromamap,0)) 
Multithread transx=trans-transy*dimsize(chromamap,0) 
string sChromatoD=name+"_14chroma"
make/O/D/N=(dimsize(trans,0),3) $sChromatoD=0
wave ChromatoD=$sChromatoD
Multithread ChromatoD[][1]=mass[transx[x]] 
Multithread ChromatoD[][2]=transz[x]
// Décalage temporel due au offset
variable timeoffset=dimoffset(chromamap,1)
findlevel/Q chromatime, timeoffset
variable timeoffsetindx=ceil(V_levelX)
Multithread ChromatoD[][0]=chromatime[transy[x]+timeoffsetindx]
NewGizmo/N=$name
AppendToGizmo axes=defaultAxes,name=axes0
ModifyGizmo displayLastObject
AppendToGizmo path=root:$sChromatoD,name=path0
ModifyGizmo displayLastObject
ModifyGizmo stopUpdates
ModifyGizmo ModifyObject=path0,objectType=path,property={ pathColorType,3}
ModifyGizmo ModifyObject=path0,objectType=path,property={ pathCTab,ColdWarm}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,gridType,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,gridPlaneColor,0,1.5259e-05,0.2,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabel,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabel,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelText,"Retention time"}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelText,"m/z"}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelCenter,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelCenter,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelDistance,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelDistance,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelScale,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelScale,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelRGBA,0.000000,0.000000,0.000000,1.000000}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelRGBA,0.000000,0.000000,0.000000,1.000000}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelTilt,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelTilt,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelFont,"default"}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelFont,"default"}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,axisLabelFlip,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,axisLabelFlip,0}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 0,labelBillboarding,1}
ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 1,labelBillboarding,1}
ModifyGizmo resumeUpdates
if(stats==1)
	// Extraction des indicateurs
	transx=ChromatoD[x][0] //tR
	transy=ChromatoD[x][1] // m/z
	transz=ChromatoD[x][2] // I
	wavestats/Q transz
	print "tR : "+num2str(transx[V_maxRowLoc])
	print "m/z : "+num2str(transy[V_maxRowLoc])
	print "Imax : "+num2str(wavemax(transz))
endif
end

////// CHROMAFIELD - WINDOWS AND HOOK FUNCTIONS

function InitChromaMapBuilder(name)
string name
string graphname="ChromaMapBuilder_"+name
if(wintype(graphname))
	killwindow $graphname
endif
string sTic=name+"_2chroma"
string stempus=name+"_3chroma"
string sMasse=name+"_4chroma"
string sSpecIndex=name+"_7chroma"
string sCropChroma=name+"_8chroma"
string sCropMassIndx=name+"_9chroma"
string sFitResults=name+"_11chroma"
string sLeanmap=name+"_12chroma"
string smapgenparam=name+"_13chroma"
wave Masse=$sMasse
wave specindex=$sSpecIndex
wave tic=$stic
wave tempus=$stempus
wave Leanmap=$sLeanmap
wave maxtimerange
wave mapgenparam=$smapgenparam
if(waveexists(mapgenparam))
	maxtimerange[0]=mapgenparam[0]
	maxtimerange[1]=mapgenparam[1]
else
	maxtimerange[0]=wavemin(tempus)
	maxtimerange[1]=wavemax(tempus)
endif
make/O/D/N=(1) NoiseLevelAjust
if(waveexists(mapgenparam))
	NoiseLevelAjust=mapgenparam[4]
else
	NoiseLevelAjust=3
endif
make/O/D/N=(1) DisplayRes4Ionmap
DisplayRes4Ionmap=0
if(waveexists(MouseHookEnabler)==0)
	make/O/D/N=(1) MouseHookEnabler
	MouseHookEnabler=0
else
	wave MouseHookEnabler
endif
make/O/D/N=(1) Int2Parents
Int2Parents=0
make/O/D/N=(1) Lean2IonMap
Lean2IonMap=0
Display/K=1/N=$graphname
ModifyGraph/W=$graphname margin(right)=60
controlbar/W=$graphname 70
Appendtograph/W=$graphname/L=TimeTic/B=Tic/Vert tic vs tempus
ModifyGraph/W=$graphname axisEnab(Tic)={0.9,1},axisEnab(TimeTic)={0,0.85},freePos(Tic)=0,freepos(TimeTic)=0,nticks(Tic)=2,fStyle(TimeTic)=1,lblPos(TimeTic)=80,lblLatPos(TimeTic)=0
Label/W=$graphname TimeTic "Retention Time"
setaxis/W=$graphname TimeTic, wavemin(tempus), wavemax(tempus)
setwindow $graphname, hook(hName)=ChromaMapBuilderHook, hookcursor=2
doupdate/W=$graphname
SetControlsChromaMap(name)
if(waveexists(Leanmap))
	mapgenerator(name,"mass")
	string sthemap=name+"_5chroma"
	wave themap=$sthemap
	appendimage/W=$graphname/L=timetic/B=MassIndex themap
	ModifyGraph/W=$graphname axisEnab(MassIndex)={0,0.89}
	ModifyImage/W=$graphname $sthemap ctab= {*,wavemax(themap)/100,ColdWarm,0}
	MatrixOP/O $sCropChroma=sumcols(themap)^t
	wave CropChroma=$sCropChroma
	findlevel/Q tempus, dimoffset(themap,1)
	variable indxbas=V_levelX
	findlevel/Q tempus, dimoffset(themap,1)+dimdelta(themap,1)*dimsize(themap,1)
	variable indxhaut=V_levelX
	variable normal=(wavemax(tic,indxbas,indxhaut)-wavemin(tic,indxbas,indxhaut))/(wavemax(CropChroma)-wavemin(CropChroma))
	variable baseline=wavemin(CropChroma)
	CropChroma=CropChroma-baseline
	CropChroma=(CropChroma)*normal
	baseline=wavemin(tic,indxbas,indxhaut)
	CropChroma+=baseline
	setscale/P x, dimoffset(themap,1), dimdelta(themap,1), CropChroma
	Appendtograph/W=$graphname/L=TimeTic/B=Tic/Vert $sCropChroma
	ModifyGraph/W=$graphname mode($sCropChroma)=7,hbFill($sCropChroma)=4,zColor($sCropChroma)={$sCropChroma,*,*,ColdWarm,0}
	variable MaxIndx, MinIndx
	MaxIndx=dimsize(themap,0)
	MinIndx=0
	make/O/D/N=(4) poslabel
		poslabel[0]=floor(((MaxIndx-MinIndx)/5)*1+MinIndx)
		poslabel[1]=floor(((MaxIndx-MinIndx)/5)*2+MinIndx)
		poslabel[2]=floor(((MaxIndx-MinIndx)/5)*3+MinIndx)
		poslabel[3]=floor(((MaxIndx-MinIndx)/5)*4+MinIndx)
	make/O/T/N=(4) loclabel
		loclabel[0]=num2str(Masse[poslabel[0]])
		loclabel[1]=num2str(Masse[poslabel[1]])
		loclabel[2]=num2str(Masse[poslabel[2]])
		loclabel[3]=num2str(Masse[poslabel[3]])
	islandfit4ionmap(name)
endif
if(waveexists(specindex))
	appendtograph/W=$graphname/L=Spec/B=MassIndex specindex
	ModifyGraph/W=$graphname axisEnab(Spec)={0.85,1}
	Label/W=$graphname Spec "Index\rSpectrum"
	ModifyGraph/W=$graphname tkLblRot(Spec)=90,freePos(Spec)=3,fStyle(Spec)=1,nticks(Spec)=2
	Label/W=$graphname MassIndex "Mass/Charge"
	ModifyGraph/W=$graphname userticks(MassIndex)={poslabel,loclabel}
	ModifyGraph/W=$graphname freePos(MassIndex)=0,lblPos(MassIndex)=40
	MatrixOP/O $sCropMassIndx=sumrows(themap)
	wave CropMassIndx=$sCropMassIndx
	duplicate/O CropMassIndx wave1stripped,dataindex,roi1
	duplicate/O Masse wave0stripped,roi0, wave2stripped, roi2
	wave2stripped[]=defmass(wave0stripped[x])
	roi2[]=defmass(wave0stripped[x])
	Makeindex/R wave1stripped, dataindex
	calcinfo()
endif
if(waveexists($sFitResults))
	Movewindow/W=$graphname 530.25,42.5,1089,539.75
	InitFitManPanel(name)
else
	Movewindow/W=$graphname 830.25,42.5,1389,539.75
endif
end

function SetControlsChromaMap(name)
string name
wave Sliderint2parent, NoiseLevelAjust
wave MouseHookEnabler,DisplayRes4Ionmap,Int2Parents,Lean2IonMap
string graphname="ChromaMapBuilder_"+name
string leanmap=name+"_12chroma"
getwindow $graphname wsizeDC
variable width=V_right-V_left, height=V_bottom-V_top
CorrectCsrChromaMap(name)
SetVariable ChromaBuilderNoise title="Set noise level for chromatogram treatment: ",pos={width/2+5,5},size={width/2-5,40},proc=NoiseLevelAjustVariable, value=NoiseLevelAjust[0],win=$graphname//, limits={0,15,0.5}
if(waveexists($leanmap))
	Button ChromaBuilderComplete title="Retreat chromatogram experiment", pos={0,0},size={width/2,70}, proc=ButtonChromaBuilderCompleteProc,win=$graphname
	CheckBox ChromaCheckMouseHook proc=ChromaBoxMouseHook,value=MouseHookEnabler[0],win=$graphname,value=MouseHookEnabler[0],pos={3*width/4-125,30}
	CheckBox ChromaCheckFitResults4Ionmap proc=ChromaBoxFitResults4Ionmap,win=$graphname,value=DisplayRes4Ionmap[0],pos={3*width/4-125,50}
	CheckBox ChromaCheckInt2Parents proc=ChromaBoxInt2Parents,win=$graphname,value=Int2Parents[0],pos={3*width/4,30}
	CheckBox ChromaCheckLean2IonMap proc=ChromaBoxLean2IonMap,win=$graphname,value=Lean2IonMap[0],pos={3*width/4,50}
	if(MouseHookEnabler[0]==0)
		CheckBox ChromaCheckMouseHook title="Activate Mouse hook",win=$graphname
	else
		CheckBox ChromaCheckMouseHook title="Desactivate Mouse hook",win=$graphname
	endif
	if(DisplayRes4Ionmap[0]==0)
		CheckBox ChromaCheckFitResults4Ionmap title="Hide fit results",win=$graphname
	else
		CheckBox ChromaCheckFitResults4Ionmap title="Display fit results",win=$graphname
	endif
	if(Int2Parents[0]==0)
		CheckBox ChromaCheckInt2Parents title="Display Island map",win=$graphname
	else
		CheckBox ChromaCheckInt2Parents title="Display Intensity map",win=$graphname
	endif
	if(Lean2IonMap[0]==0)
		CheckBox ChromaCheckLean2IonMap title="Display Ion Map",win=$graphname
	else
		CheckBox ChromaCheckLean2IonMap title="Display Lean Map",win=$graphname
	endif
	TitleBox CursorValue title="Cursor value",pos={width-95+30,75},frame=0,win=$graphname
	ValDisplay CursortRDisplay title="tR",pos={width-85,90},size={80,20},bodyWidth=50,value=_NUM:Nan,win=$graphname
	ValDisplay CursormzDisplay title="m/z",pos={width-85,105},size={80,20},bodyWidth=50,value=_NUM:Nan,win=$graphname
	ValDisplay CursorIDDisplay title="ID",pos={width-85,120},size={80,20},bodyWidth=50,value=_NUM:Nan,win=$graphname
else
	Button ChromaBuilderComplete title="Treat chromatogram experiment", pos={0,0},size={width/2,70}, proc=ButtonChromaBuilderCompleteProc,win=$graphname
endif
end

Function ChromaMapBuilderHook(s)
	STRUCT WMWinHookStruct &s
	string graphname=s.winname
	string name=ReplaceString("ChromaMapBuilder_",s.winname,"")
	string chromaname="ChromatoViewPanel_"+name
	string panelname="FitManPanel_"+name
	string sthemap=name+"_5chroma"
	string sMasse=name+"_4chroma"
	string sTic=name+"_2chroma"
	string sSpecIndex=name+"_7chroma"
	string stempus=name+"_3chroma"
	string sCropChroma=name+"_8chroma"
	string sCropMassIndx=name+"_9chroma"
	string sParentMap=name+"_10chroma"
	string sFitResults=name+"_11chroma"
	string sLeanmap=name+"_12chroma"
	string smapgenparam=name+"_13chroma"
	string schromatoD=name+"_14chroma"
	//string sHKmap=name+"_10chroma"
	string sActivemap=Replacestring(";",ImageNameList( s.winname,";"),"")
	wave/Z CropChroma=$sCropChroma
	wave/Z CropMassIndx=$sCropMassIndx
	wave/Z tempus=$stempus
	wave/Z specindex=$sSpecIndex
	wave/Z tic=$stic
	wave/Z masse=$smasse
	wave/Z activemap=$sActivemap
	wave/Z leanmap=$sLeanmap
	wave/Z mapgenparam=$smapgenparam
	if(waveexists($sthemap))
		wave themap=$sthemap
		wave wave1stripped,wave0stripped,wave2stripped,roi1,roi2,roi0
		variable mouseIndex=axisvalfrompixel(s.winname, "MassIndex", s.mouseloc.h)
		variable mouseTime=axisvalfrompixel(s.winname, "TimeTic", s.mouseloc.v)
		variable mouseID=Nan
		variable mouseMasse=Masse[mouseIndex]
		GetAxis/W=$s.winname/Q TimeTic
		variable inTime=((mouseTime<=V_max)&&(mouseTime>=V_min))
		GetAxis/W=$s.winname/Q MassIndex
		variable inIndx=((mouseindex<=V_max)&&(mouseindex>=V_min))
		getmarquee/W=$s.winname
	endif
	
	switch( s.eventCode )
		case 2: //Window killed
			// Place for cleaning stuff
			if(wintype(chromaname))
				killwindow $chromaname
			endif
			if(wintype(panelname))
				killwindow $panelname
			endif
			KillFitList(name)
			removeimage/Z/W=$graphname $sParentMap
			removeimage/Z/W=$graphname $sthemap
			RemoveFromGraph/Z/W=$graphname fitres4ionmap
			RemoveFromGraph/Z/W=$graphname $sCropChroma
			RemoveFromGraph/Z/W=$graphname $sCropMassIndx
			killwaves/Z $sthemap,$sParentMap,fitres4ionmap,$sCropChroma,$sCropMassIndx,$schromatoD
			killwaves/Z islandchroma,pli,extrindx,extrtime,extrpli,extr,persi,death,birth,fit_islandchroma,point2move,trace2move,poslabel,loclabel // ADD CW 13032019
			break
		case 3: //mouse down
			break
		case 4: //mouse moved
			// A refaire - problèmes hors des axes (intempestif)
			wave MouseHookEnabler
			if(MouseHookEnabler[0]==1)
				if(waveexists($sthemap))
					Cursor/A=0/H=0/S=0/C=(0,0,0)/W=$s.winname C $sSpecIndex mouseIndex
					findlevel/Q tempus, mouseTime
					Cursor/A=0/H=0/S=0/C=(0,0,0)/P/W=$s.winname D $sTic V_levelX
					GetAxis/W=$s.winname/Q Spec
					ValDisplay IndxBox title="m/z=",pos={s.mouseloc.h,PixelFromAxisVal(s.winname, "Spec",V_max)+56},value=_NUM:mouseMasse,bodyWidth=60,frame=2, win=$s.winname,mode=3, disable=0
					GetAxis/W=$s.winname/Q Tic
					ValDisplay TimeBox title="RT=",pos={pixelfromaxisval(s.winname,"Tic",V_max),s.mouseloc.v+43},value=_NUM:mouseTime,bodyWidth=40,frame=2, win=$s.winname,mode=3, disable=0
					ValDisplay IntBox title="Intensity=",pos={0,70},value=_NUM:themap(mouseindex)(mousetime),bodyWidth=80,frame=2, win=$s.winname,mode=3, disable=0
				endif
			endif
			break
		case 5: //mouse up
			if(waveexists($sthemap))
				variable Csrtr,Csrmz,CsrID
				if(inTime==0 && inIndx==1) //Hors temps et en masse
					variable currTime=max(dimoffset(themap,1)+dimdelta(themap,1)*str2num(stringbykey("YPOINT",csrinfo(E,s.winname))),dimoffset(themap,1))
					Cursor/I/H=1/S=2/C=(0,52224,52224)/W=$s.winname E $sActivemap mouseIndex,currTime
					CropChroma=themap(mouseIndex)(x)
					ModifyGraph/W=$graphname mode($sCropChroma)=7,hbFill($sCropChroma)=4,zColor($sCropChroma)={$sCropChroma,*,wavemax(themap)/100,ColdWarm,1}
					Duplicate/O SpecIndex wave1stripped
					Duplicate/O Masse wave0stripped, wave2stripped
					wave2stripped[]=defmass(wave0stripped[x])
					Csrtr=currTime
					Csrmz=mouseMasse
					// Extraction des colonnes
					extract/O/FREE/INDX leanmap,transindx,q==0
					extract/O/FREE leanmap,transmass,q==0
					extract/O/FREE leanmap,transtime,q==1
					// Extraction des masses correspondantes
					extract/O/FREE transindx, transindx,transmass==floor(mouseIndex)
					extract/O/FREE transtime, transtime,transmass==floor(mouseIndex)
					// Extraction des temps correspondants
					extract/O/FREE transindx, transindx,transtime==floor(currtime/dimdelta(themap,1)-dimoffset(themap,1)/dimdelta(themap,1))
					CsrID=leanmap[transindx[0]][2]*(numtype(transindx[0])==0)
				elseif(inTime==1 && inIndx==0) // En temps et hors masse
					variable currIndx=str2num(stringbykey("POINT",csrinfo(E,s.winname)))
					Cursor/I/H=1/S=2/C=(0,52224,52224)/W=$s.winname E $sActivemap currIndx,mouseTime
					CropMassIndx=themap(x)(mousetime)
					Duplicate/O CropMassIndx roi1
					Duplicate/O Masse roi0
					Csrtr=max(min(mouseTime,mapgenparam[1]),mapgenparam[0])
					if(currIndx<0 || numtype(currindx)==2)
						Csrmz=Nan
						CsrID=Nan
					else
						Csrmz=Masse[currIndx]
						// Extraction des colonnes
						extract/O/FREE/INDX leanmap,transindx,q==0
						extract/O/FREE leanmap,transmass,q==0
						extract/O/FREE leanmap,transtime,q==1
						// Extraction des masses correspondantes
						extract/O/FREE transindx, transindx,transmass==floor(currIndx)
						extract/O/FREE transtime, transtime,transmass==floor(currIndx)
						// Extraction des temps correspondants
						extract/O/FREE transindx, transindx,transtime==floor(Csrtr/dimdelta(themap,1)-dimoffset(themap,1)/dimdelta(themap,1))
						CsrID=leanmap[transindx[0]][2]*(numtype(transindx[0])==0)
					endif 
				elseif(inTime==1 && inIndx==1) //En temps et en masse
					Cursor/I/H=1/S=2/C=(0,52224,52224)/W=$s.winname E $sActivemap mouseIndex,mouseTime
					CropMassIndx=themap(x)(mousetime)
					CropChroma=themap(mouseIndex)(x)
					ModifyGraph/W=$graphname mode($sCropChroma)=7,hbFill($sCropChroma)=4,zColor($sCropChroma)={$sCropChroma,*,wavemax(themap)/100,ColdWarm,1}
					Duplicate/O CropMassIndx roi1
					Duplicate/O Masse roi0
					Csrtr=max(min(mouseTime,mapgenparam[1]),mapgenparam[0])
					Csrmz=mouseMasse
					// Extraction des colonnes
					extract/O/FREE/INDX leanmap,transindx,q==0
					extract/O/FREE leanmap,transmass,q==0
					extract/O/FREE leanmap,transtime,q==1
					// Extraction des masses correspondantes
					extract/O/FREE transindx, transindx,transmass==floor(mouseIndex)
					extract/O/FREE transtime, transtime,transmass==floor(mouseIndex)
					// Extraction des temps correspondants
					extract/O/FREE transindx, transindx,transtime==floor(Csrtr/dimdelta(themap,1)-dimoffset(themap,1)/dimdelta(themap,1))
					CsrID=leanmap[transindx[0]][2]*(numtype(transindx[0])==0)
				else // Hors temps et hors masse => Reinitialisation du curseur & affichages
					Cursor/K/W=$s.winname E
					MatrixOP/O $sCropChroma=sumcols(themap)^t
					wave CropChroma=$sCropChroma
					findlevel/Q tempus, dimoffset(themap,1)
					variable indxbas=V_levelX
					findlevel/Q tempus, dimoffset(themap,1)+dimdelta(themap,1)*dimsize(themap,1)
					variable indxhaut=V_levelX
					variable normal=(wavemax(tic,indxbas,indxhaut)-wavemin(tic,indxbas,indxhaut))/(wavemax(CropChroma)-wavemin(CropChroma))
					variable baseline=wavemin(CropChroma)
					CropChroma=CropChroma-baseline
					CropChroma=(CropChroma)*normal
					baseline=wavemin(tic,indxbas,indxhaut)
					CropChroma+=baseline
					setscale/P x, dimoffset(themap,1), dimdelta(themap,1), CropChroma
					ModifyGraph/W=$graphname mode($sCropChroma)=7,hbFill($sCropChroma)=4,zColor($sCropChroma)={$sCropChroma,*,*,ColdWarm,0}
					Csrtr=Nan
					Csrmz=Nan
					CsrID=Nan
				endif
				// Affichage des valeurs
				ValDisplay CursortRDisplay title="tR", value=_NUM:(Csrtr),win=$graphname
				ValDisplay CursormzDisplay title="m/z", value=_NUM:(Csrmz),win=$graphname
				ValDisplay CursorIDDisplay title="ID", value=_NUM:(CsrID),win=$graphname
				// Affichage des options par clic droit
				if(s.eventMod==16 || s.eventMod==17) //(Left click) or (Left click + right click)
					ChromatoContextualMenu(CsrID,name)
				endif
				// Surlignage des lignes de fit dans la fenêtre appropriée
				string Fitlist="fitlist_"+name
				string Fitsel="selmode_"+name
				wave/T listwave=$Fitlist
				wave selmode=$Fitsel
				wave fitresults=$sfitresults
				resortDATAbycol(selmode,listwave,1)
				make/O/D/N=(dimsize(fitresults,0)) transID
				transID=fitresults[x][0]
				string transnum=listwave[0][1]
				if(str2num(transnum)!=wavemin(transID))
					resortDATAbycol(selmode,listwave,1)
				endif
				extract/O/D/INDX/FREE fitresults,transfitresults,transID==CsrID
				variable k=0,n=dimsize(transfitresults,0)
				selmode[][0]=64 //RAZ des sélections
				for(k=0;k<n;k+=1)
					selmode[transfitresults[k]][0]+=1
				endfor
				ListBox fitlist row=-1,win=$panelname
			endif
			Getaxis/W=$graphname/Q MassIndex
			wave poslabel
				poslabel[0]=floor(((V_max-V_min)/5)*1+V_min)
				poslabel[1]=floor(((V_max-V_min)/5)*2+V_min)
				poslabel[2]=floor(((V_max-V_min)/5)*3+V_min)
				poslabel[3]=floor(((V_max-V_min)/5)*4+V_min)
			wave/T loclabel
				loclabel[0]=num2str(Masse[poslabel[0]])
				loclabel[1]=num2str(Masse[poslabel[1]])
				loclabel[2]=num2str(Masse[poslabel[2]])
				loclabel[3]=num2str(Masse[poslabel[3]])										   
			break
		case 6: //Resize
			SetControlsChromaMap(name)
			Getaxis/W=$graphname/Q MassIndex
			wave poslabel
				poslabel[0]=floor(((V_max-V_min)/5)*1+V_min)
				poslabel[1]=floor(((V_max-V_min)/5)*2+V_min)
				poslabel[2]=floor(((V_max-V_min)/5)*3+V_min)
				poslabel[3]=floor(((V_max-V_min)/5)*4+V_min)
			wave/T loclabel
				loclabel[0]=num2str(Masse[poslabel[0]])
				loclabel[1]=num2str(Masse[poslabel[1]])
				loclabel[2]=num2str(Masse[poslabel[2]])
				loclabel[3]=num2str(Masse[poslabel[3]])
			break
		case 7: //cursor moved
			string csrname=s.cursorname
			wave maxtimerange
				if(stringmatch(csrname,"A"))
					maxtimerange[0]=hcsr($csrname)
				elseif(stringmatch(csrname,"B"))
					maxtimerange[1]=hcsr($csrname)
				endif
			CorrectCsrChromaMap(name)
			break
		case 12://moved
			// Exemple : AdvMolManPanel qui bouge avec ses fenêtres graphiques et stoechio
			//if(wintype(graphname))
				//autopositionwindow/M=0/R=$panelname $graphname
			//endif
			break
		case 22: //mouse wheel
			break
	endswitch
	return 0
end

function ChromatoContextualMenu(CsrID,name)
variable CsrID
string name
string sfitresults=name+"_11chroma"
wave fitresults=$sfitresults
// Génération de _10chroma
string hkmap=name+"_10chroma"
mapgenerator(name,"island")
PopupContextualMenu "Automatic - Initial guesses;Automatic - Fit this island;Manual - Change guesses;Manual - New guesses;Manual - Fit this island;Save this island fit;Kill this island fit"
strswitch(S_selection)	
	case "Automatic - Initial guesses":
		if(CsrID!=Nan)
			island2chroma(CsrID,name)
			pershomo(islandchroma)
			Fitguessfrompershomo(persi,extrtime,0)
			fitislandguess(islandchroma, fitinit)
		endif
		break
	case "Automatic - Fit this island":
		if(CsrID!=Nan)
			island2chroma(CsrID,name)
			pershomo(islandchroma)
			Fitguessfrompershomo(persi,extrtime,0)
			fitisland(islandchroma, fitinit)
		endif
		break
	case "Manual - Change guesses":
		if(CsrID!=Nan)
			island2chroma(CsrID,name)
			pershomo(islandchroma)
		endif
		break
	case "Manual - New guesses":
		if(CsrID!=Nan)
			island2chroma(CsrID,name)
			Fitguessfrompershomo(persi,extrtime,1)
			fitislandguess(islandchroma, fitinit)
		endif
		break
	case "Manual - Fit this island":
		if(CsrID!=Nan)
			island2chroma(CsrID,name)
			Fitguessfrompershomo(persi,extrtime,1)
			fitisland(islandchroma, fitinit)
		endif
		break
	case "Save this island fit":
		// Suppresion de l'ancien fit pour l'ilot
		make/O/D/FREE/N=(dimsize(fitresults,0)) transID
		transID=fitresults[x][0]
		extract/O/D/INDX/FREE fitresults,transfitresults,transID==CsrID
		variable meanmass=fitresults[transfitresults[0]][1]
		deletepoints/M=0 transfitresults[0],numpnts(transfitresults),fitresults
		// Ajout du nouveau fit pour l'ilot
		wave fitcoef, fitinit
		variable numfit=numpnts(fitcoef)/4,l,modelcount=0
		insertpoints/M=0 0,numfit,fitresults
		for(l=0;l<numfit*4;l+=4)
			fitresults[modelcount][0]=CsrID //IlotID
			fitresults[modelcount][1]=meanmass //mean mass of the island
			// Paramètres après fit
			fitresults[modelcount][2]= fitcoef[l]//H
			fitresults[modelcount][3]= abs(fitcoef[l+1])//c
			fitresults[modelcount][4]= abs(fitcoef[l+2])//l
			fitresults[modelcount][5]= fitcoef[l+3]//exp
			// Valeurs d'initialisation
			fitresults[modelcount][6]= fitinit[l]//H
			fitresults[modelcount][7]= fitinit[l+1]//c
			fitresults[modelcount][8]= abs(fitinit[l+2])//l
			fitresults[modelcount][9]= fitinit[l+3]//exp
			// Calcul de l'aire
			fitresults[modelcount][10]=fitarea4island(name,islandchroma,modelcount)
			// Passage à la ligne suivante
			modelcount+=1
		endfor
		// Tri selon les IslandID de _11chroma
		variable col=0
		variable n=dimsize(fitresults,0)
		Make/FREE/O/D/N=(n) key
		Make/FREE/O/D/N=(n) index
		key=fitresults[x][col]+x
		Makeindex/A key, index
		Duplicate/FREE/O fitresults newfitresults
		if(index[0]<index[numpnts(key)-1])
   			wavetransform/O flip index
   			if(fitresults[index[0]][col]>fitresults[index[n-1]][col])
   				wavetransform/O flip index
   			endif
		endif
		newfitresults=fitresults[index[x]][y]
		Duplicate/O newfitresults fitresults
		// MAJ de la liste
		MajFitList(name)
		break
	case "Kill this island fit":
		// Suppresion de l'ancien fit pour l'ilot
		make/O/D/N=(dimsize(fitresults,0)) transID
		transID=fitresults[x][0]
		extract/O/D/INDX fitresults,transfitresults,transID==CsrID
		deletepoints/M=0 transfitresults[0],numpnts(transfitresults),fitresults
		MajFitList(name)
		break
endswitch
killwaves $hkmap
end

function InitFitManPanel(name)
string name
InitFitList(name)
string mapname="ChromaMapBuilder_"+name
string chromaname="ChromatoViewPanel_"+name
string panelname="FitManPanel_"+name
string Fitlist="fitlist_"+name
string Fitsel="selmode_"+name
string Fittitle="fittitle_"+name
string mapgens=name+"_13chroma"
wave mapgen=$mapgens
if(wintype(chromaname))
	killwindow $chromaname
endif
if(wintype(panelname))
	killwindow $panelname
endif
MAKE/O/N=(0) islandchroma,pli,extrindx,extrtime,extrpli,extr,persi,death,birth,fit_islandchroma,point2move // ADD CW 13032019
MAKE/O/T/N=(0) trace2move // ADD CW 14032019
//Afficher ici chromato, tic et fit
Display/K=1/N=$chromaname 
ModifyGraph/W=$chromaname margin(right)=20,margin(left)=30,margin(top)=20,margin(bottom)=40
Movewindow/W=$chromaname 0,42.5,200,539.75
autopositionwindow/M=0/R=$mapname $chromaname
controlbar/W=$chromaname 70
Appendtograph/W=$chromaname/L=TimeTic/B=Tic/Vert islandchroma
Appendtograph/W=$chromaname/L=TimeTic/B=Tic/Vert persi vs extrtime
Appendtograph/W=$chromaname/L=TimeTic/B=Tic/Vert fit_islandchroma
ModifyGraph/W=$chromaname rgb(islandchroma)=(0,0,0),mode(persi)=8,marker(persi)=19, rgb(fit_islandchroma)=(0,0,65535)
GetAxis/W=$chromaname/Q TimeTic
ModifyGraph/W=$chromaname freePos(TimeTic)={0,Tic}, freePos(Tic)={V_min,TimeTic}
SetAxis/W=$chromaname/A/N=1 TimeTic
ModifyGraph/W=$chromaname lblPos(Tic)=35,lblPos(TimeTic)=45,standoff=0,axisEnab(Tic)={0.07,1},axisEnab(TimeTic)={0.025,0.865}
Label/W=$chromaname Tic "Intensity"
Label/W=$chromaname TimeTic "Retention Time"
getwindow $chromaname wsizeDC
variable width=V_right-V_left
Button ChromatoViewMarkerButton title="Remove markers", pos={0,0},size={width/2,70/2}, proc=ButtonChromatoViewRemoveFitMarkersProc,win=$chromaname
Button ChromatoViewShowGuessButton title="Show guess", pos={0,35},size={width/2,70/2}, proc=ButtonChromatoViewShowGuessProc,win=$chromaname
Button ChromatoViewNewFitButton title="Process fit", pos={width/2,0},size={width/2,70/2}, proc=ButtonChromatoViewNewFitProc,win=$chromaname
Button ChromatoViewSaveFitButton title="Save this fit", pos={width/2,35},size={width/2,70/2}, proc=ButtonChromatoViewSaveFitProc,win=$chromaname
ValDisplay ChromatoViewIDDisplay title="Island ID:",pos={40,75},size={100,30},bodyWidth=80,value=_NUM:Nan,win=$chromaname
ValDisplay ChromatoViewMZDisplay title="Peak mass:",pos={40,95},size={100,30},bodyWidth=80,value=_NUM:Nan,win=$chromaname
ValDisplay ChromatoViewNoiseDisplay title="Noise cut:",pos={40,115},size={100,30},bodyWidth=80,value=_NUM:mapgen[5],win=$chromaname
SetWindow $chromaname hook(dragname)=ChromatoViewPanelHook
//Afficher ici liste des fits
Newpanel/K=1/N=$panelname/W=(0,42.5,517,539.75) 
autopositionwindow/M=0/R=$chromaname $panelname
getwindow $panelname wsizeDC
width=V_right-V_left
variable height=V_bottom-V_top
ListBox fitlist,pos={1,1},size={width-15,height-15},proc=PilotFitList,win=$panelname
ListBox fitlist listWave=$Fitlist,selWave=$Fitsel,titleWave=$Fittitle,mode=4,widths={8,20,20,20,20,20,20,20,20},userColumnResize=1,win=$panelname
SetWindow $panelname hook(hookname)=FitManPanelHook
end

function ChromatoViewPanelHook(s)
	STRUCT WMWinHookStruct &s
	string name=ReplaceString("ChromatoViewPanel_",s.winname,"")
	string chromaname="ChromatoViewPanel_"+name
	wave persi,extrtime,point2move
	wave/T trace2move
	variable res=0
	variable lepoint=str2num(stringbykey("HITPOINT",TraceFromPixel(s.mouseloc.h,s.mouseloc.v,"")))
	variable xx,yy
	switch(s.eventCode)
		case 3:// mouse down
			res=1
			if(s.eventmod==1) // Left clic : move a point
				insertpoints 0,1,point2move,trace2move
				point2move[0]=lepoint
				trace2move[0]="persi"
			elseif(s.eventmod==9)
				make/O/N=(0) persi,extrtime
				//print "Ctrl+click"
			elseif(s.eventmod==17) // Right clic : add a point
				insertpoints 0,1,point2move,trace2move,persi,extrtime
				point2move[0]=0
				trace2move[0]="persi"
				//ADD CW 20092019 - Add the point at the right click location
				xx=axisvalfrompixel(chromaname,"Tic",s.mouseloc.h)
				yy=axisvalfrompixel(chromaname,"TimeTic",s.mouseloc.v)
				getAxis/W=$chromaname/Q Tic
				persi[point2move[0]]=limit(xx,V_min,V_max)
				extrtime[point2move[0]]=yy
			endif
			break
		case 4:// mouse moved
		if(strlen(trace2move[0]))
			xx=axisvalfrompixel(chromaname,"Tic",s.mouseloc.h)
			yy=axisvalfrompixel(chromaname,"TimeTic",s.mouseloc.v)
			GetAxis/W=$chromaname/Q TimeTic
			if(yy<V_min || yy >V_max) // Delete point : go outiside the time limits
				deletepoints point2move[0],1,extrtime,persi
				trace2move="" //RAZ to avoid erase of all points
			else // Move point
				getAxis/W=$chromaname/Q Tic
				persi[point2move[0]]=limit(xx,V_min,V_max)
				extrtime[point2move[0]]=yy
			endif
			res=0
		endif
			break
		case 5: //mouse up
		deletepoints 0,1, point2move,trace2move
		sort extrtime,extrtime,persi
			break
	endswitch
	return res
end

Function FitManPanelHook(s)
	STRUCT WMWinHookStruct &s
	string panelname=s.winname
	switch( s.eventCode )
		case 6: //Resize
			//controlinfo/W=panelname 
			getwindow $panelname wsizeDC
			variable width=V_right-V_left, height=V_bottom-V_top
			ListBox fitlist size={width-15,height-15}
			break
	endswitch
	return 0
end

Function PilotFitList(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
	string name=ReplaceString("FitManPanel_",lba.win,"")
	string graphname="ChromaMapBuilder_"+name
	string chromaname="ChromatoViewPanel_"+name
	string hkmap=name+"_10chroma"
	string sfitresults=name+"_11chroma"
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
	wave time2indx=$name + "_3chroma"//temps <-> time index
	wave mass2indx=$name + "_4chroma" //mass <-> mass index
	wave intmap=$name + "_5chroma"//int map
	wave hkresults=$name + "_6chroma"//HK table
	wave leanmap=$name + "_12chroma"//lean map
	wave mapgenparam=$name + "_13chroma"//map generation parameters
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 2:
			if(row==-1) // column sorting
				resortDATAbycol(selwave,listwave,col)
			endif
			if(col==0)//Contextual actions
				ChromatoContextualMenu(str2num(lba.listWave[row][1]),name)
			endif
			break
		case 3: // double click
			// Génération de _10chroma					
			mapgenerator(name,"island")
			// Génération du chromato et positions de pics automatique
			island2chroma(str2num(lba.listWave[row][1]),name)
			pershomo(islandchroma)
			Fitguessfrompershomo(persi,extrtime,0)
			// Reprise des paramètres initiaux de fit enregistrés dans _11chroma					  
			wave fitresults=$sfitresults
			make/O/D/N=(dimsize(fitresults,0)) transID
			transID=fitresults[x][0]
			extract/O/D/INDX transID,transfitresults,transID==str2num(lba.listWave[row][1])
			variable k=0,fitcount=0,n=dimsize(transfitresults,0)*4
			make/O/D/N=(n) fitinit
			for(k=0;k<n;k+=4)
				//H
				fitinit[k]=fitresults[transfitresults[fitcount]][6]
				//c
				fitinit[k+1]=fitresults[transfitresults[fitcount]][7]
				//l
				fitinit[k+2]=fitresults[transfitresults[fitcount]][8]
				//exp
				fitinit[k+3]=1
				// Compteur de fit
				fitcount+=1
			endfor
			// Fit
			fitisland(islandchroma, fitinit)
			// Interactions with chromatoview
			variable islandid4view, mass4view
			islandid4view=str2num(lba.listWave[row][1])
			mass4view=str2num(lba.listWave[row][2])
			GetAxis/Q/W=$graphname TimeTic
			SetAxis/W=$chromaname TimeTic V_min,V_max
			ValDisplay ChromatoViewIDDisplay value=_NUM:(islandid4view),win=$chromaname
			ValDisplay ChromatoViewMZDisplay value=_NUM:(mass4view),win=$chromaname
			killwaves $hkmap
			break
		case 4: // cell selection
			// Zoom autour de l'ilot sur la carte
			variable islandnum=dimsize(hkresults,0)
			make/O/D/FREE/N=(islandnum) transnum
			transnum=hkresults[x][8]
			// on trouve la ligne dans HK results de l'ilot en question
			findlevel/Q transnum, str2num(lba.listWave[row][1]) //**//
			if(V_flag==1)
				Abort("Error - Island not found\r\rIsland not found error could be due to:\r- Dirac in mass;\r- Less than 10 points in time;\r- Too large mass range;\r- The background.")
			else 
				variable ilotline=round(V_levelX)
			endif
			// Offset de la map en temps
			variable minmass = hkresults[ilotline][2]
			variable maxmass = hkresults[ilotline][3]
			variable mintime = hkresults[ilotline][4]
			variable maxtime = hkresults[ilotline][5]
			minmass=findmassortimevalue(minmass,mass2indx)//*min mass indx*//
			maxmass=findmassortimevalue(maxmass,mass2indx)//*max mass indx*//
			SetAxis/W=$graphname MassIndex minmass-(maxmass-minmass)*0.1,maxmass+(maxmass-minmass)*0.1
			SetAxis/W=$graphname TimeTic mintime-(maxtime-mintime)*0.1,maxtime+(maxtime-mintime)*0.1
			SetControlsChromaMap(name)
			string sMasse=name+"_4chroma"
			wave/Z masse=$smasse
			Getaxis/W=$graphname/Q MassIndex
			wave poslabel
				poslabel[0]=floor(((V_max-V_min)/5)*1+V_min)
				poslabel[1]=floor(((V_max-V_min)/5)*2+V_min)
				poslabel[2]=floor(((V_max-V_min)/5)*3+V_min)
				poslabel[3]=floor(((V_max-V_min)/5)*4+V_min)
			wave/T loclabel
				loclabel[0]=num2str(Masse[poslabel[0]])
				loclabel[1]=num2str(Masse[poslabel[1]])
				loclabel[2]=num2str(Masse[poslabel[2]])
				loclabel[3]=num2str(Masse[poslabel[3]])
			//expectamaximavectorbuilder(hkresults[ilotline][8], name)
			break
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 12://keystroke
			//print row
			if(row==127)
				wave fitresults=$sfitresults
				// Suppresion de l'ancien fit pour les cellules sélectionnées
				extract/O/INDX/FREE selwave,killme,selwave[x][y]>64 && y==0
				if(dimsize(killme,0)==0)
					variable ID2kill=str2num(listWave[killme[0]][1])
					make/O/D/FREE/N=(dimsize(fitresults,0)) transID
					transID=fitresults[x][0]
					extract/O/D/INDX/FREE fitresults,transfitresults,transID==ID2kill
					deletepoints/M=0 transfitresults[0],numpnts(transfitresults),fitresults
				else
					variable j, num2kill=dimsize(killme,0)
					for(j=0;j<num2kill;j+=1)
						ID2kill=str2num(listWave[killme[j]][1])
						make/O/D/FREE/N=(dimsize(fitresults,0)) transID
						transID=fitresults[x][0]
						extract/O/D/INDX/FREE fitresults,transfitresults,transID==ID2kill
						if(dimsize(transfitresults,0)!=0)
							deletepoints/M=0 transfitresults[0],numpnts(transfitresults),fitresults
						endif
					endfor
				endif
				MajFitList(name)
				resortDATAbycol(selwave,listwave,2)
			elseif(row==8)
				// Génération de _10chroma
				mapgenerator(name,"island")
				// Génération du chromato et positions de pics automatique
				extract/O/INDX/FREE selwave,killme,selwave[x][y]>64 && y==0
				variable ID2graph=str2num(listWave[killme[0]][1])
				variable mass2graph=str2num(listWave[killme[0]][2])
				island2chroma(ID2graph,name)
				pershomo(islandchroma)
				Fitguessfrompershomo(persi,extrtime,0)
				// Reprise des paramètres initiaux de fit enregistrés dans _11chroma
				wave fitresults=$sfitresults
				make/O/D/N=(dimsize(fitresults,0)) transID
				transID=fitresults[x][0]
				extract/O/D/INDX transID,transfitresults,transID==ID2graph
				k=0
				fitcount=0
				n=dimsize(transfitresults,0)*4
				make/O/D/N=(n) fitinit
				for(k=0;k<n;k+=4)
					//H
					fitinit[k]=fitresults[transfitresults[fitcount]][6]
					//c
					fitinit[k+1]=fitresults[transfitresults[fitcount]][7]
					//l
					fitinit[k+2]=fitresults[transfitresults[fitcount]][8]
					//exp
					fitinit[k+3]=1
					// Compteur de fit
					fitcount+=1
				endfor
				// Fit
				fitisland(islandchroma, fitinit)
				//print sqrt(variance(islandchroma))/mean(islandchroma)*100
				// Interactions with chromatoview
				GetAxis/Q/W=$graphname TimeTic
				SetAxis/W=$chromaname TimeTic V_min,V_max
				ValDisplay ChromatoViewIDDisplay value=_NUM:(ID2graph),win=$chromaname
				ValDisplay ChromatoViewMZDisplay value=_NUM:(mass2graph),win=$chromaname
				killwaves $hkmap
			endif
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch
	return 0
End

function expectamaximavectorbuilder(ilotID, samplename)
variable ilotID
string samplename
wave time2indx=$samplename + "_3chroma"//temps <-> time index
wave mass2indx=$samplename + "_4chroma" //mass <-> mass index
wave ionmap=$samplename + "_5chroma" //ion map
wave hkresults=$samplename + "_6chroma"//HK table
// Génération de _10chroma					
mapgenerator(samplename,"island")
wave hkmap=$samplename + "_10chroma" //HK map
wave mapgenparam=$samplename + "_13chroma"//map generation parameters

variable islandnum=dimsize(hkresults,0)
make/O/D/FREE/N=(islandnum) transnum
multithread transnum=hkresults[x][8]
// on trouve la ligne dans HK results de l'ilot en question
findlevel/Q transnum, ilotID //**//
if(V_flag==1)
	Abort("Error - Island not found\r\rIsland not found error could be due to:\r- Dirac in mass;\r- Less than 10 points in time;\r- Too large mass range;\r- The background.")
else 
	variable ilotline=round(V_levelX)
endif

wave meanmass
// Offset de la map en temps
variable mapoffset=mapgenparam[0] //min time of the map
variable minmass = hkresults[ilotline][2]
variable maxmass = hkresults[ilotline][3]
variable mintime = hkresults[ilotline][4]-mapoffset
variable maxtime = hkresults[ilotline][5]-mapoffset
print minmass, maxmass,mintime,maxtime

minmass=findmassortimevalue(minmass,mass2indx)//*min mass indx*//
maxmass=findmassortimevalue(maxmass,mass2indx)//*max mass indx*//
mintime=findmassortimevalue(mintime,time2indx)//*min time indx*//
maxtime=findmassortimevalue(maxtime,time2indx)//*max time indx*//
// Islandchroma = vecteur 1D qui somme les intensités pour différentes masses appartenant à l'ilotID
make/O/D/N=((maxtime-mintime)*(maxmass-minmass)) Xcol,Ycol
variable l,m,counter=0
for(l=minmass;l<=maxmass;l+=1)
	for(m=mintime;m<=maxtime;m+=1)
		Xcol[counter]=(l)*(hkmap[l][m]==ilotID)
		Ycol[counter]=(m)*(hkmap[l][m]==ilotID)
		counter+=1
		//islandchroma[m-mintime]+=ionmap[l][m]*(hkmap[l][m]==ilotID)
	endfor
endfor
extract/O Xcol, Xcol, Xcol!=0
extract/O Ycol, Ycol, Ycol!=0
Xcol=mass2indx[Xcol[x]]
Ycol=time2indx[Ycol[x]]
killwaves hkmap
end

Function InitFitList(name)
string name
string sselmode="selmode_"+name
string sfitlist="fitlist_"+name
string sfittitle="fittitle_"+name
string sfitresults=name+"_11chroma"
wave fitresults=$sfitresults
variable n=dimsize(fitresults,0)
Make/O/I/N=(n,7) $sselmode=0
Make/O/T/N=(n,7) $sfitlist=""
Make/O/T/N=7 $sfittitle=""
wave selmode=$sselmode
wave/T fitlist=$sfitlist
wave/T fittitle=$sfittitle
// Intitulés de colonne
fittitle={" ","Island ID","Peak Mass","Peak RT","Peak height","Peak Area","Peak Shape"}
// Type de cellule : 64=triangle de selection, 32=carré à cocher
selmode[][0]=64 
// Construction de la liste
variable k=0
string s1,s2,s3,s4,s5,s6,s7,s8,s9
for(k=0;k<n;k+=1)
	sprintf s1,"%.0f",fitresults[k][0] //IslandID
	fitlist[k][1]=s1
	sprintf s2,"%.3f",fitresults[k][1] //Peak mass
	fitlist[k][2]=s2
	sprintf s3,"%.2f",fitresults[k][3] //Peak tR
	fitlist[k][3]=s3
	sprintf s4,"%.0f",fitresults[k][2] // Peak height
	fitlist[k][4]=s4
	sprintf s5,"%.0f",fitresults[k][10] // Peak area
	fitlist[k][5]=s5
	sprintf s6,"%.2f",1/(fitresults[k][4]*fitresults[k][5]) // Peak asymetry
	fitlist[k][6]=s6
endfor
end

Function MajFitList(name)
string name
string sselmode="selmode_"+name
string sfitlist="fitlist_"+name
string sfitresults=name+"_11chroma"
wave fitresults=$sfitresults
variable n=dimsize(fitresults,0)
Make/O/I/N=(n,7) $sselmode=0
Make/O/T/N=(n,7) $sfitlist=""
wave selmode=$sselmode
wave/T fitlist=$sfitlist
// Type de cellule : 64=triangle de selection, 32=carré à cocher
selmode[][0]=64 
// Construction de la liste
variable k=0
string s1,s2,s3,s4,s5,s6,s7,s8,s9
for(k=0;k<n;k+=1)
	sprintf s1,"%.0f",fitresults[k][0] //IslandID
	fitlist[k][1]=s1
	sprintf s2,"%.3f",fitresults[k][1] //Peak mass
	fitlist[k][2]=s2
	sprintf s3,"%.2f",fitresults[k][3] //Peak tR
	fitlist[k][3]=s3
	sprintf s4,"%.0f",fitresults[k][2] // Peak height
	fitlist[k][4]=s4
	sprintf s5,"%.0f",fitresults[k][10] // Peak area
	fitlist[k][5]=s5
	sprintf s6,"%.2f",1/(fitresults[k][4]*fitresults[k][5]) // Peak asymetry
	fitlist[k][6]=s6
endfor
listelesdonnees()
end

function KillFitList(name)
string name
string selmode="selmode_"+name
string fitlist="fitlist_"+name
string fittitle="fittitle_"+name
killwaves/Z $selmode,$fitlist,$fittitle,fitres4ionmap
end

function CorrectCsrChromaMap(name)
string name
string graphname="ChromaMapBuilder_"+name
string sTic=name+"_2chroma"
string stempus=name+"_3chroma" 
string smapgenparam=name+"_13chroma"
wave mapgenparam=$smapgenparam
wave maxtimerange
sort maxtimerange,maxtimerange
if(maxtimerange[0]<wavemin($stempus))
	maxtimerange[0]=wavemin($stempus)
endif
if(maxtimerange[1]>wavemax($stempus))
	maxtimerange[1]=wavemax($stempus)
endif
GetAxis/W=$graphname/Q TimeTic
variable Csrmin=V_min,  Csrmax=V_max
GetAxis/W=$graphname/Q Tic
if((maxtimerange[0]>=Csrmin)&&(maxtimerange[1]<=Csrmax)) // MODIF CW 19032019
	cursor/F/N=1/C=(34816,34816,34816)/H=3/S=2 A $stic maxtimerange[0],wavemax($stempus)
	cursor/F/N=1/C=(34816,34816,34816)/H=3/S=2 B $stic maxtimerange[1],wavemax($stempus)
	SetVariable A title=" ",pos={pixelfromaxisval(graphname,"Tic",V_max),63+pixelfromaxisval(graphname,"Timetic",maxtimerange[0])},value=maxtimerange[0],proc=SetTimeRangeA,limits={wavemin($stempus),wavemax($stempus),1},win=$graphname,disable=0
	SetVariable B title=" ",pos={pixelfromaxisval(graphname,"Tic",V_max),63+pixelfromaxisval(graphname,"Timetic",maxtimerange[1])},value=maxtimerange[1],proc=SetTimeRangeB,limits={wavemin($stempus),wavemax($stempus),1},win=$graphname,disable=0
elseif((maxtimerange[0]<Csrmin)&&(maxtimerange[1]<=Csrmax)) // MODIF CW 19032019
	cursor/K A
	cursor/F/N=1/C=(34816,34816,34816)/H=3/S=2 B $stic maxtimerange[1],wavemax($stempus)
	SetVariable A win=$graphname,disable=1
	SetVariable B title=" ",pos={pixelfromaxisval(graphname,"Tic",V_max),63+pixelfromaxisval(graphname,"Timetic",maxtimerange[1])},value=maxtimerange[1],proc=SetTimeRangeB,limits={wavemin($stempus),wavemax($stempus),1},win=$graphname,disable=0
elseif((maxtimerange[0]>=Csrmin)&&(maxtimerange[1]>Csrmax)) // MODIF CW 19032019
	cursor/F/N=1/C=(34816,34816,34816)/H=3/S=2 A $stic maxtimerange[0],wavemax($stempus)
	cursor/K B
	SetVariable A title=" ",pos={pixelfromaxisval(graphname,"Tic",V_max),63+pixelfromaxisval(graphname,"Timetic",maxtimerange[0])},value=maxtimerange[0],proc=SetTimeRangeA,limits={wavemin($stempus),wavemax($stempus),1},win=$graphname,disable=0
	SetVariable B win=$graphname,disable=1
else
	cursor/K A
	cursor/K B
	SetVariable A win=$graphname,disable=1
	SetVariable B win=$graphname,disable=1
endif
end


Function SetTimeRangeA(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
		string name=replacestring("ChromaMapBuilder_",sva.win,"")
			Variable dval = sva.dval
			String sval = sva.sval
			string stempus=name+"_3chroma"
			wave tempus=$stempus
			Findvalue/T=(dimdelta(tempus,0)/2)/V=(dval) tempus
			wave maxtimerange
			maxtimerange[0]=tempus[V_value]
			CorrectCsrChromaMap(name)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetTimeRangeB(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
		string name=replacestring("ChromaMapBuilder_",sva.win,"")
			Variable dval = sva.dval
			String sval = sva.sval
			string stempus=name+"_3chroma"
			wave tempus=$stempus
			Findvalue/T=(dimdelta(tempus,0)/2)/V=(dval) tempus
			wave maxtimerange
			maxtimerange[1]=tempus[V_value]
			CorrectCsrChromaMap(name)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonChromaBuilderCompleteProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
string name=ReplaceString("ChromaMapBuilder_",ba.win,"")
	switch( ba.eventCode )
		case 2: // mouse up
			setwindow $ba.win, hook(hname)=$""
			print "****************"
			print "Chromafield module"
			print "Starting data treatment"
			print "Sample name : ",name
			ConstruitMap(name)
			ChromatoMap2vectorsHK(name)
			//MainEMchroma(name)
			island2chromastruct(name)
			InitChromaMapBuilder(name)
			listelesdonnees()
			print "****************"
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ChromaBoxLean2IonMap(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	wave Lean2IonMap
	string name=ReplaceString("ChromaMapBuilder_",cba.win,"")
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Lean2IonMap=checked
			if(checked==0)
				mapgenerator(name,"mass")
				CheckBox ChromaCheckLean2IonMap title="Display Ion Map"
			else
				// Récupération des paramètres de génération précédents de la carte
				wave maxtimerange,NoiseLevelAjust
				string smapgenparam=name+"_13chroma"
				wave mapgenparam=$smapgenparam
				maxtimerange[0]=mapgenparam[0]
				maxtimerange[1]=mapgenparam[1]
				NoiseLevelAjust=mapgenparam[4]
				//Construction de la carte
				ConstruitMap(name)
				CheckBox ChromaCheckLean2IonMap title="Display Lean Map"
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ChromaBoxMouseHook(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	wave MouseHookEnabler
	string name=ReplaceString("ChromaMapBuilder_",cba.win,"")
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			MouseHookEnabler=checked
			InitChromaMapBuilder(name)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ChromaBoxFitResults4Ionmap(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	wave DisplayRes4Ionmap
	string name=ReplaceString("ChromaMapBuilder_",cba.win,"")
	string graphname="ChromaMapBuilder_"+name
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			DisplayRes4Ionmap=checked
			if(DisplayRes4Ionmap[0]==1)
				RemoveFromGraph/W=$graphname/Z fitres4ionmap
				CheckBox ChromaCheckFitResults4Ionmap title="Display fit results"
			else
				islandfit4ionmap(name)
				CheckBox ChromaCheckFitResults4Ionmap title="Hide fit results"
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

// Permet de supperposer des points sur l'ion map
function islandfit4ionmap(name)
string name
string graphname="ChromaMapBuilder_"+name
wave mass2indx=$name + "_4chroma" //mass <-> mass index
wave islandres=$name + "_11chroma" //modelisation results
make/O/D/N=(dimsize(islandres,0),3) fitres4ionmap
fitres4ionmap[][0]=islandres[x][1]
fitres4ionmap[][1]=islandres[x][3]
fitres4ionmap[][2]=islandres[x][7]
variable mass,k,n=dimsize(islandres,0)
for(k=0;k<n;k+=1)
	mass = fitres4ionmap[k][0]
	findlevel/Q mass2indx, mass //*mass indx*//
	if(V_flag==1)
		mass=0
	else 
		mass=round(V_levelX)
	endif
	fitres4ionmap[k][0]=mass
endfor
appendtograph/W=$graphname/L=timetic/B=massindex fitres4ionmap[][1] vs fitres4ionmap[][0]
ModifyGraph/W=$graphname mode(fitres4ionmap)=3,marker(fitres4ionmap)=19,rgb(fitres4ionmap)=(65535,65535,65535)
end

Function NoiseLevelAjustVariable(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	wave NoiseLevelAjust
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			NoiseLevelAjust[0]=str2num(sval)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonChromatoViewRemoveFitMarkersProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			make/O/N=(0) persi,extrtime
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function ButtonChromatoViewShowGuessProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			print "TBD"
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function ButtonChromatoViewNewFitProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			print "TBD"
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function ButtonChromatoViewSaveFitProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			print "TBD"
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

function Fitted2spectrum(name)
string name
string sFit=name+"_11chroma"
wave Fitted=$sFit
string sint="roi1"
string smass="roi0"
wave int=$sint
wave mass=$smass
extract/O Fitted,int,y==2 //&& Fitted!=0
extract/O Fitted,mass,y==1 //&& Fitted!=0
Duplicate/O mass roi2
wave roi0=$smass
roi2=defmass(roi0)
end

Function ButtonShowFittedPeaks(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
		controlinfo/W=AdvancedManager TabCon_List_CHROMA
		wave/t lasel=$S_value
		string name=lasel[V_value][1]+"_11chroma"
		Fitted2spectrum(lasel[V_value][1])
		//edit/K=1 $name
			break
	endswitch
	return 0
End

Function ButtonDuplicateChroma(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=advancedmanager TabCon_List_CHROMA
			// Name of the active chroma
			wave/T lasel=$S_value
			string name=lasel[V_value][1]
			// Ask for a new name
			string newname
			prompt newname, "New Name"
			Doprompt "Enter a new good name for this sample",newname
			if(v_flag==0)
				if(cmpstr(name,newname)==0)
					print "The two names are not different"
				else
					renameactivechroma(name, newname)
				endif
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function renameactivechroma(name, newname)
string name, newname
string s0chroma=name+"_0chroma"
string s1chroma=name+"_1chroma"
string s2chroma=name+"_2chroma"
string s3chroma=name+"_3chroma"
string sn0chroma=newname+"_0chroma"
string sn1chroma=newname+"_1chroma"
string sn2chroma=newname+"_2chroma"
string sn3chroma=newname+"_3chroma"
Duplicate/O $s0chroma, $sn0chroma
Duplicate/O $s1chroma, $sn1chroma
Duplicate/O $s2chroma, $sn2chroma
Duplicate/O $s3chroma, $sn3chroma
listelesdonnees()
end

Function ChromaBoxInt2Parents(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	string graphname=cba.win
	string name=ReplaceString("ChromaMapBuilder_",cba.win,"")
	string sleanmap=name+"_12chroma"
	string sthemap=name+"_5chroma"
	string sParentMap=name+"_10chroma"
	wave leanmap=$sleanmap
	wave Int2Parents
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Int2Parents=checked
			if(checked==0)
				mapgenerator(name,"mass")
				wave/Z themap=$sthemap
				wave/Z parentmap=$sParentMap
				removeimage/Z/W=$graphname $sParentMap
				removeimage/Z/W=$graphname $sthemap
				appendimage/W=$graphname/L=TimeTic/B=MassIndex themap
				ModifyGraph/W=$graphname axisEnab(MassIndex)={0,0.89}
				ModifyImage/W=$graphname $sthemap ctab= {*,wavemax(themap)/100,ColdWarm,0}
				killwaves/Z parentmap
				CheckBox ChromaCheckInt2Parents title="Display Island map"
			else
				mapgenerator(name,"island")
				wave/Z themap=$sthemap
				wave/Z parentmap=$sParentMap
				removeimage/Z/W=$graphname $sParentMap
				removeimage/Z/W=$graphname $sthemap
				appendimage/W=$graphname/L=TimeTic/B=MassIndex parentmap
				ModifyGraph/W=$graphname axisEnab(MassIndex)={0,0.89}
				ModifyImage/W=$graphname $sParentMap ctab={*,*,Classification,0}
				killwaves/Z themap
				CheckBox ChromaCheckInt2Parents title="Display Intensity map"
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function Findtr4mz(mz, ppmprecision, themap, timetic, specmass, specint,smoothmode, windowsearch)
	variable mz, ppmprecision, smoothmode, windowsearch
	wave themap, timetic, specmass, specint
	variable adjprecision=ppmprecision*mz*10^(-6)
	// On cherche l'index en masse le plus proche de notre mz renseigné
	FindValue/V=(mz)/T=(adjprecision) specmass
	if(V_Value==-1)
		print "No mass at the specified value"
		return -1
	endif
	// On identifie le max d'intensité en masse autour de ce point
	variable maxint=wavemax(specint, pnt2x(specint,V_Value-windowsearch),pnt2x(specint,V_Value+windowsearch))
	findlevel/Q/R=(pnt2x(specint,V_Value-windowsearch),pnt2x(specint,V_Value+windowsearch)) specint, maxint
	// On identifie le max d'intensité en temps relatif au max d'intensité en masse
	Extract/O/FREE themap, transchroma, p==V_LevelX
	if(smoothmode==1)
		Smooth/B 13, transchroma
	endif
	variable ciblemass=specmass[V_LevelX]
	variable retime=wavemax(transchroma)
	// Décalage temporel due au offset
	variable timeoffset=dimoffset(themap,1)
	findlevel/Q timetic, timeoffset
	variable timeoffsetindx=ceil(V_levelX)
	// On extrait le tR
	FindValue/V=(retime) transchroma
	retime=timetic[V_Value+timeoffsetindx]
	//retime=timetic[V_Value]
	print "[m/z;tR]=["+num2str(ciblemass)+";"+num2str(retime)+"]"
end

////// END CHROMAFIELD - WINDOWS AND HOOK FUNCTIONS

// 12/11/2018 CW: fonction to kill null intensities in spectra waves
function killzeroint(swavename)
string swavename
wave zerostripped=$swavename+"_0stripped"
wave onestripped=$swavename+"_1stripped"
// Index des intensités non nulles
extract/O/FREE/INDX onestripped, transindx, onestripped!=0
// On refait les waves
make/O/FREE/D/N=(dimsize(transindx,0)) transx,transy
transx=zerostripped[transindx[x]]
transy=onestripped[transindx[x]]
duplicate/O transx, zerostripped
duplicate/O transy, onestripped
end

//////////////////////////////////////////////////////////////
///////// END CHROMATOGRAMS
//////////////////////////////////////////////////////////////

Function Mendev_Formula_Set(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
wave molecule
nvar charge
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			//Duplicate/O str2mol(sval) molecule
			molecule=str2mol(sval)[x]
			charge = str2charge(sval)
			tuelesspec()
			genestringmol()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

ThreadSafe Function target2theo(target,error)
variable target,error
return (-1e6*target)/(error-1e6)
end

Function MolbagRescale(name)
string name
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
wave wave1stripped,wave0stripped
variable k,n=dimsize(ciblebag,0),masse
print n
for(k=0;k<n;k+=1)
	//extract/O/FREE molbag,molecule, y==k
	masse=target2theo(ciblebag[k][1],ciblebag[k][3])//mol2massmono(molecule,chargebag[k])//
	ciblebag[k][1]=wave0stripped[prochepic(ciblebag[k][1],inf)]
	ciblebag[k][2]=wave1stripped[prochepic(ciblebag[k][1],inf)]
	ciblebag[k][3]=1e6*(masse-ciblebag[k][1])/masse
endfor
end

function MolbagRename(name,newname)
string name,newname
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
smolbag="Molbag_"+newname
schargebag="Chargebag_"+newname
slistbag="Listbag_"+newname
sciblebag="Ciblebag_"+newname
Rename molbag,$smolbag
Rename chargebag,$schargebag
Rename listbag,$slistbag
Rename ciblebag,$sciblebag
end

///////////////////////////////////////////////////////////
//////////GRAPHTTRIBUTOR////////////////////
///////////////////////////////////////////////////////////

/////INTERFACE////////

Function initGRAPHTTRIBUTOR(name)
string name
string graphname="GRAPHTTRIBUTOR"
if(wintype(graphname))
	killwindow $graphname
endif
Display/K=1/N=$graphname
Movewindow/W=$graphname 20,20,200,300
controlbar 50
setwindow $graphname, hook(hName )=GRAPHTTRIBUTORHook, hookcursor=2
make/O/D/N=(2,0) Vertice2BagSettings
Vertice2BagSettings[0]=10
Vertice2BagSettings[1]=5
variable/G verticemin,seedppmthreshold
verticemin=Vertice2BagSettings[0]
seedppmthreshold=Vertice2BagSettings[1]
InitVariator(name)
InitResGraph(0)
SetControlsGRAPHTTRIBUTOR(name)
formatStds()
end

Function InitVariator(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=dimsize(molbag,1)
make/O/D/N=(n,3) VariatorSel
make/O/T/N=(n,3) VariatorList
VariatorSel=48*(y==0)
VariatorList[][0]=""
VariatorList[][1]=listbag[x]
VariatorList[][2]=num2str(target2theo(ciblebag[x][1],ciblebag[x][3]))
end

Function InitResGraph(nb)
variable nb
Make/O/T/N=(3) ResGraphTitle
ResGraphTitle={"Min. mass node","Vertices","Max. Int."}
Make/O/D/N=(nb,3) ResGraphSel=0
Make/O/T/N=(nb,3) ResGraphList="test"
end

function SetControlsGRAPHTTRIBUTOR(name)
string name
string graphname="GRAPHTTRIBUTOR"
wave VariatorSel,ResGraphSel
wave/T VariatorList,ResGraphTitle,ResGraphList
variable n=dimsize(VariatorSel,0)
getwindow $graphname wsizeDC
variable width=V_right-V_left, height=V_bottom-V_top, nb=3
PopupMenu popup1, win=$graphname,pos={0,49},bodyWidth=width,proc=PopMenuGraphVari,mode=30,popvalue=(name),value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
ListBox mollist,win=$graphname,pos={0,71},size={width,max(50,n*24)},labelBack=(47872,47872,47872),frame=0,listWave=VariatorList,selWave=VariatorSel,widths={15,75,75},userColumnResize= 0,special= {0,20,0}, proc=VariatorSelListAction
ListBox treelist,win=$graphname,pos={0,71+max(50,n*24)},size={width,height-21-max(50,n*24)},mode=1,frame=0,titlewave=ResGraphTitle,listWave=ResGraphList,selWave=ResGraphSel,widths={25,25,25},userColumnResize= 0,proc=TreeListAction
Button Forest, win=$graphname,pos={0,0}, size={width/nb,50}, proc=ButtonForest,title="Grow Forest"
Button Vertice2Bag, win=$graphname, pos={0+width/nb,0}, size={width/nb,50}, proc=ButtonVertice2Bag, title="Forest to Bag"
end

Function ButtonForest(ba)
	STRUCT WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			GrowForest()
			break
	endswitch
	return 0
end

Function ButtonVertice2Bag(ba)
	STRUCT WMButtonAction &ba
	string lesvertice2bag=winlist("Vertice2BagPanel*",";",""), levertice2bag
	variable k, n=itemsinlist(lesvertice2bag)
	switch( ba.eventCode )
		case 2: // mouse up
			if(stringmatch(levertice2bag,""))
				Execute "Vertice2BagPanel()"
			else
				for(k=0;k<n;k+=1)
					levertice2bag=stringfromlist(k,lesvertice2bag)
					killwindow $levertice2bag
				endfor
			endif
			break
	endswitch
end

Window Vertice2BagPanel()
	PauseUpdate; Silent 1		// building window...
	NewPanel /K=1/W=(188,425,418,530)
	autopositionwindow/M=0/R=GRAPHTTRIBUTOR Vertice2BagPanel
	SetVariable setvar1,pos={10,10},size={200,16},title="Consider all vertices down to..."
	SetVariable setvar1,value=verticemin,live= 1
	SetVariable setvar2,pos={10,32},size={200,16},title="With this seed threshold (ppm)"
	SetVariable setvar2,value=seedppmthreshold,live= 1
	Button button1,pos={64,60},size={100,40},proc=GoForVertice2Bag,title="Go!"
	Button button1,fColor=(65280,0,0)
end

Function GoForVertice2Bag(ba)
	STRUCT WMButtonAction &ba
	nvar verticemin,seedppmthreshold
	switch(ba.eventcode)
		case 2:
			string bagname="Vertice2Bag"
			prompt bagname, "Enter a name for the molbag destination"
			Doprompt "Enter a list name for the newly attributed items to land in",bagname
			AttribSelTribue(verticemin,seedppmthreshold,bagname)
			break
	endswitch
	return 0
end

// Function ADDED CW 16042019 for generation of a molbag from selected graphtribution list
function AttribSelTribue(numvertices, tolerance,bag)
variable numvertices, tolerance
string bag
// Extraction index of the tribue with more or equal to numvertices
wave/T ResGraphList
wave Ciblebag_Propositions
wave Adja_Obs,wave0stripped,wave1stripped,Stds_Obs_Mol2add,Stds_Obs_charge2add
variable k, vertice
vertice=str2num(ResGraphList[0][1])
for(k=0;vertice>=numvertices;k+=1)
	PropagAttrib(Adja_Obs,wave0stripped,wave1stripped,str2num(ResGraphList[k][0]),"Propagation",Stds_Obs_Mol2add,Stds_Obs_charge2add)
	if(abs(Ciblebag_Propositions[0][3])<=tolerance)
		MolbagDecant("Propagation",bag,Molbag2FullIndx("Propagation"))
	endif
	//print "test"+num2str(k)+"& "+num2str(vertice)
	vertice=str2num(ResGraphList[k+1][1])
endfor
end

function AttribSelTribueAsym(numvertices, tolerancemax, tolerancemin, bag)
variable numvertices, tolerancemax, tolerancemin
string bag
// Extraction index of the tribue with more or equal to numvertices
wave/T ResGraphList
wave Ciblebag_Propositions
wave Adja_Obs,wave0stripped,wave1stripped,Stds_Obs_Mol2add,Stds_Obs_charge2add
variable k, vertice
vertice=str2num(ResGraphList[0][1])
for(k=0;vertice>=numvertices;k+=1)
	PropagAttrib(Adja_Obs,wave0stripped,wave1stripped,str2num(ResGraphList[k][0]),"Propagation",Stds_Obs_Mol2add,Stds_Obs_charge2add)
	if(Ciblebag_Propositions[0][3]<=tolerancemax && Ciblebag_Propositions[0][3]>=tolerancemin)
		MolbagDecant("Propagation",bag,Molbag2FullIndx("Propagation"))
	endif
	//print "test"+num2str(k)+"& "+num2str(vertice)
	vertice=str2num(ResGraphList[k+1][1])
endfor
end

Function TreeListAction(la)
	STRUCT WMListboxAction &la
	switch( la.eventCode )
		case 2: // mouse up
			if(la.row==-1)
				resortDATAbycol(la.selwave,la.listwave,la.col)
			else
				wave Parent_Obs
				Arbre2ROI(Parent_Obs,str2num(la.listWave[la.row][0]))
				//speClu4chromato(0)
				//fast3Droi()
				PropagAttrib(Adja_Obs,wave0stripped,wave1stripped,str2num(la.listWave[la.row][0]),"Propagation",Stds_Obs_Mol2add,Stds_Obs_charge2add)
			endif
		case -1: // control being killed
			break
	endswitch
	return 0
end

Function VariatorSelListAction(la)
	STRUCT WMListboxAction &la
	switch( la.eventCode )
		case 2: // mouse up
			formatStds()
		case -1: // control being killed
			break
	endswitch
	return 0
end

Function PopMenuGraphVari(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			InitVariator(popStr)
			SetControlsGRAPHTTRIBUTOR(popStr)
			formatStds()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function GRAPHTTRIBUTORHook(s)
	STRUCT WMWinHookStruct &s
	ControlInfo /W=$s.winName popup1
	switch( s.eventCode )
		case 2:
			CleanUpGRAPHTTRIBUTOR()
			break
		case 3:
			break
		case 4:
			break
		case 5:
			break
		case 6:
			 SetControlsGRAPHTTRIBUTOR(S_value)
			 break
		case 7:
			break
	endswitch
	return 0
end

Function CleanUpGRAPHTTRIBUTOR()
string graphname="GRAPHTTRIBUTOR"
Make/o/d/n=0 Traceur_obs
KillWaves/Z VariatorSel, VariatorList,ResGraphTitle,ResGraphSel,ResGraphList,Adja_obs,Poids_obs,Traceur_obs,Degree_Obs,Parent_Obs,Chef_Obs,Tribue_Obs,Stds_Obs_Mol2add,Stds_Obs,Std_Obs_Charge2add
end

function storeGraph(name)
string name
string dest
wave ResGraphSel,ResGraphList,Adja_obs,Poids_obs,Traceur_obs,Degree_Obs,Parent_Obs,Chef_Obs,Tribue_Obs,Stds_Obs_Mol2add,Stds_Obs,Std_Obs_Charge2add
dest="ResGraphSel_"+name
duplicate/O ResGraphSel, $dest
dest="ResGraphList_"+name
duplicate/O ResGraphList, $dest
dest="Adja_"+name
duplicate/O Adja_obs, $dest
dest="Poids_"+name
duplicate/O Poids_obs, $dest
dest="Traceur_"+name
duplicate/O Traceur_obs, $dest
dest="Degree_"+name
duplicate/O Degree_Obs, $dest
dest="Parent_"+name
duplicate/O Parent_Obs, $dest
dest="Chef_"+name
duplicate/O Chef_Obs, $dest
dest="Tribue_"+name
duplicate/O Tribue_Obs, $dest
dest="Stds_Mol2add_"+name
duplicate/O Stds_Obs_Mol2add, $dest
dest="Stds_"+name
duplicate/O Stds_Obs, $dest
dest="Std_Charge2add_"+name
duplicate/O Std_Obs_Charge2add, $dest
end

function loadGraph(name)
string name
string dest
wave ResGraphSel,ResGraphList,Adja_obs,Poids_obs,Traceur_obs,Degree_Obs,Parent_Obs,Chef_Obs,Tribue_Obs,Stds_Obs_Mol2add,Stds_Obs,Std_Obs_Charge2add
dest="ResGraphSel_"+name
duplicate/O $dest, ResGraphSel
dest="ResGraphList_"+name
duplicate/O $dest, ResGraphList
dest="Adja_"+name
duplicate/O $dest,Adja_obs
dest="Poids_"+name
duplicate/O $dest, Poids_obs
dest="Traceur_"+name
duplicate/O $dest, Traceur_obs
dest="Degree_"+name
duplicate/O $dest, Degree_Obs
dest="Parent_"+name
duplicate/O $dest, Parent_Obs
dest="Chef_"+name
duplicate/O $dest, Chef_Obs
dest="Tribue_"+name
duplicate/O $dest, Tribue_Obs
dest="Stds_Mol2add_"+name
duplicate/O $dest, Stds_Obs_Mol2add
dest="Stds_"+name
duplicate/O $dest, Stds_Obs
dest="Std_Charge2add_"+name
duplicate/O $dest, Std_Obs_Charge2add
end

////////METIER///////

function/WAVE initiateAdja(Obs)
wave obs
variable n=dimsize(obs,0)
make/O/D/N=(n,n)/FREE Adja
//Adja=abs(x-y)==1
return Adja
end

Function/WAVE calcBestof(Obs,Adja,Stds)//need sorted Obs et 1-dim, modifie Adja
wave obs,Adja,stds
variable n=dimsize(obs,0)
variable dim=dimsize(obs,1)
variable trivmax=wavemax(Stds)+1, ecart
make/O/D/N=(n,n)/FREE distances=Inf
Duplicate/O/FREE Stds, Candidates
variable k,j,contdim
for(k=0;k<n;k+=1)
	for(j=k+1;j<n;j+=1)
		ecart=abs(obs[k]-obs[j])
		if(ecart<trivmax)
			Candidates=abs(ecart-stds)
			if(wavemin(Candidates)<0.5)
				extract/indx/free candidates, whichgap, candidates==wavemin(Candidates)
				distances[k][j]=wavemin(Candidates)
				distances[j][k]=wavemin(Candidates)
				Adja[k][j]=whichgap[0]+1
				Adja[j][k]=whichgap[0]+1
			else
				Adja[k][j]=0
				Adja[j][k]=0
			endif
		else
			Adja[k][j]=0
			Adja[j][k]=0
			j=n
		endif
	endfor
endfor
return Distances //sort une Adja qui est l'index du std+1
end

function/WAVE kruskal4mass(obs,Adja,Poids,Parent_Obs,degmax)//modifie Parent_Obs et retourne AdjaKruskal, degmax est le degre maximum ajouté au graph par ce kruskal
wave obs,Adja,Poids,Parent_Obs
variable degmax
duplicate/o/free Adja AdjaKruskal
AdjaKruskal=0
variable n=dimsize(obs,0)
Make/O/D/N=(n) Degree=0
extract/o poids,listepoids,x<y && numtype(poids)==0
extract/indx/o/free poids,listeindx,x<y && numtype(poids)==0
duplicate/o/free listeindx, col,row
Multithread col=floor(listeindx[x]/n)
Multithread row=listeindx[x]-col[x]*n
sort listepoids,listepoids,col,row
variable m=numpnts(listepoids),k
						   
variable stalsize=100, a
Make/FREE/O/D/N=(stalsize) stal,diff,model
stal=inverseerf(x/stalsize)
a=mean(stal)
stal/=a
duplicate/O listepoids meanpoids, residue
residue=0
meanpoids=mean(listepoids,0,x)
variable bcut=bayesianCut(listepoids)
//Make/O/D/N=(m,3) statipoids
//Statipoids[][0]=(mean(listepoids,0,x+2))/(x+2)
//Statipoids[][1]=sqrt(variance(listepoids,0,x+3))/(x+3)
//statipoids[][2]=Statipoids[x][0]-Statipoids[x][1]
m=7000
for(k=0;k<m;k+=1)//&& statipoids[k][2]>0
	model=stal*mean(meanpoids,0,k)
	diff=(model[x]-meanpoids[x*k/stalsize])
	residue[k]=(norm(diff))*k
	if(residue[k]<residue[k-1])//k==bcut)//
		print k, bcut
		//k=m
	endif
	if(FindKruskal(col[k],Parent_Obs)!=FindKruskal(row[k],Parent_Obs) && max(degree[row[k]],degree[col[k]])<degmax)
		AdjaKruskal[col[k]][row[k]]=Adja[col[k]][row[k]]
		AdjaKruskal[row[k]][col[k]]=Adja[row[k]][col[k]]
		unionKruskal(row[k],col[k],Parent_Obs)
		degree[row[k]]+=1
		degree[col[k]]+=1
	endif
endfor
//listepoids/=(x+2)
return AdjaKruskal
end

Function FindKruskal(k,parent)
Variable k
Wave parent
	If(parent[k]!=k)
		parent[k] = FindKruskal(parent[k],parent)
	endif
	return parent[k]
end

Function unionKruskal(a, b,parent)
	Variable a,b
	Wave parent
	parent[FindKruskal(a,parent)]=FindKruskal(b,parent)
End

function/WAVE buildTrace(Adja,Obs_Vizu,parent)
wave Adja,Obs_Vizu,parent
variable n=dimsize(Adja,0)
variable dim=dimsize(Obs_Vizu,1)
variable k,j
make/O/D/N=(0,dim+2)/FREE Traceur
for(k=0;k<n;k+=1)
	for(j=k+1;j<n;j+=1)
		if(Adja[j][k])
			insertpoints/m=0 0,3,traceur
			traceur[2][0,dim-1]=Obs_Vizu[k][y]
			traceur[2][dim]=findKruskal(k,Parent)
			traceur[2][dim+1]=Adja[j][k]-1
			traceur[1][0,dim-1]=Obs_Vizu[j][y]
			traceur[1][dim]=findKruskal(j,Parent)
			traceur[1][dim+1]=Adja[j][k]-1
			traceur[0][0,dim-1]=Nan
			traceur[0][dim]=Nan
			traceur[0][dim+1]=Nan
		endif
	endfor
endfor
return Traceur
end

function/WAVE buildTrace4chromato(Adja,Obs_Vizu,parent)
wave Adja,Obs_Vizu,parent
variable n=dimsize(Adja,0)
variable dim=dimsize(Obs_Vizu,1)
variable k,j
make/O/D/N=(0,dim+2)/FREE Traceur
for(k=0;k<n;k+=1)
	for(j=k+1;j<n;j+=1)
		if(Adja[j][k])
			insertpoints/m=0 0,3,traceur
			traceur[2][0,dim-1]=Obs_Vizu[k][y]
			traceur[2][dim]=findKruskal(k,Parent)
			traceur[2][dim+1]=Adja[j][k]-1
			traceur[1][0,dim-1]=Obs_Vizu[j][y]
			traceur[1][dim]=findKruskal(j,Parent)
			traceur[1][dim+1]=Adja[j][k]-1
			traceur[0][0,dim-1]=Nan
			traceur[0][dim]=Nan
			traceur[0][dim+1]=Nan
		endif
	endfor
endfor
return Traceur
end

function/Wave ClassifGridNode(Adja,Textds,Seltds,name)
wave Adja,Seltds
wave/T Textds
string name
string sname="GridNodeColor_"+name
variable n=dimsize(adja,0)//nombre de noeuds
Extract/O/FREE/T Textds,Stds,y==1 && Seltds[x][0]==48
variable m=numpnts(Stds)//nombre de liens différents
variable omega=2^m//nombre de partitions binaires
Make/O/D/N=(omega)/FREE PartiCounts=0, transcount
Make/O/T/N=(omega) PartiLabel="Empty",translabel
string inter = " ∩ "
string provname
variable k,j
for(k=1;k<omega;k+=1)
	provname=""
	for(j=0;j<m;j+=1)
		if(mod(floor(k/(2^(m-1-j))),2))
			provname=provname+inter+Stds[m-1-j]
		endif
	endfor
	PartiLabel[k]=provname[5,strlen(provname)-1]
endfor
Make/O/D/N=(m)/FREE nombre
Make/O/D/N=(n)/FREE belongsto, transbel
variable adresse
for(k=0;k<n;k+=1)
	nombre=0
	for(j=0;j<n;j+=1)
		if(adja[k][j])
			nombre[m-adja[k][j]]=2^(m-adja[k][j])
		endif
	endfor
	adresse=sum(nombre)
	Particounts[adresse]+=1
	belongsto[k]=adresse
endfor
Make/O/D/N=(omega,m)/FREE crible
crible=mod(floor(x/2^(m-1-y)),2)
MatrixOP/O/FREE Ordre=sumrows(crible)
Duplicate/O/FREE ordre, index
Makeindex ordre index
//translabel=PartiLabel[index[x]]
//transcount=Particounts[index[X]]
//transbel=index[belongsto[x]]
//belongsto=transbel[x]
//PartiCounts=transcount[x]
//PartiLabel=translabel[x]
Duplicate/O belongsto $sname
return Particounts
end

function/WAVE calcDegree(Adja)
wave Adja
matrixop/O/FREE Degree=sumrows(Adja)
return Degree
end

Function/WAVE calcChef(parent)
wave parent
Duplicate/O/FREE parent,Chef
Sort chef,chef
variable k, n=numpnts(chef)
for(k=1;k<n;k+=1)
	if(chef[k]==chef[k-1])
		deletepoints k,1,chef
		k-=1
		n-=1
	endif
endfor
Return Chef
end

function/wave CompteTribue(chef,Adja,Obs,Parent)
wave chef,adja,obs, parent
variable n=dimsize(adja,0),k
make/o/d/n=(n)/FREE passed=0
n=dimsize(chef,0)
make/o/d/n=(n,3)/FREE StatTribue=0
for(k=0;k<n;k+=1)
	Arbre2ROI(Parent,chef[k])
	wave roipnts,roi1
	StatTribue[k][0]=roipnts[0]
	StatTribue[k][1]=numpnts(roipnts)
	StatTribue[k][2]=wavemax(roi1)
endfor
return StatTribue
end

function DFScomptage(Adja,Obs,Mark,Nexus)
wave Adja, Obs, Mark
variable Nexus
variable n=dimsize(Obs,0),k
Mark[Nexus]+=1
//Action sur le sommet
//Fin d'action sur le sommet
for(k=0;k<n;k+=1)
	if(Adja[Nexus][k])
		if(Mark[k]==0)
			DFScomptage(Adja,Obs,Mark,k)
		endif
	endif
endfor
end

Function PropagAttrib(Adja,Mass,Int,Startidx,bag,std2add,charge2add)
wave Adja, mass,int,std2add,charge2add
variable Startidx
string bag
MolbagWithdrawFrom(Molbag2FullIndx(bag),bag)
nvar inco,charge
variable n=dimsize(Adja,0),k, t=ticks
inco=Mass[Startidx]
analysemass("arg")
//MolbagExtractor("Propositions",{0},bag)
InitMolList(0,"Propositions")
print "Stoichiometric initiation took ",(ticks-t)/60," seconds"
t=ticks
Duplicate/O/FREE mass, Mark
Mark=0
//Mark[Startidx]+=1
for(k=0;k<n && Adja[Startidx][k]==0;k+=1)
endfor
Extract/O/FREE molbag_propositions,molecule,y==0
DFSattrib(Adja,Mass,Int,Mark,k,molecule,Startidx,charge,bag,std2add,charge2add)
print "Stoichiometric propagation took ",(ticks-t)/60," seconds"
Molbag2generator("Propositions", 0)
end

function DFSattrib(Adja,Mass,Int,Mark,Nexus,Mothermol,motheridx,mothercharge,bag,std2add,charge2add)
wave Adja, Mass,Int, Mark, Mothermol,std2add,charge2add
variable Nexus,mothercharge,motheridx
string bag
variable n=dimsize(Adja,0),k,charge, addidx,massesimu
Mark[Nexus]+=1
//Action sur le sommet
addidx=adja[Nexus][Motheridx]-1
Duplicate/O/FREE Mothermol, molecule
molecule=Mothermol-sign(mol2massmono(mothermol,mothercharge)-Mass[Nexus])*std2add[x][addidx]
charge=mothercharge+charge2add[addidx]
massesimu=mol2massmono(molecule,charge)
MolbagAddTo(molecule,charge,mol2str(molecule,charge),{0,Mass[Nexus],int[Nexus],1e6*(massesimu-mass[nexus])/massesimu,0},bag)
//Fin d'action sur le sommet
for(k=0;k<n;k+=1)
	if(Adja[Nexus][k])
		if(Mark[k]==0)
			DFSattrib(Adja,Mass,Int,Mark,k,molecule,Nexus,charge,bag,std2add,charge2add)
		endif
	endif
endfor
end

Function Arbre2ROI(Parent,start)
wave parent
variable start
wave wave0stripped,wave1stripped,wave2stripped, wave3stripped
wave incox,incoy,incoxdm
nvar inco
inco=wave0stripped[start]
incox[0]=wave0stripped[start]
incoy[0]=wave1stripped[start]
incoxdm[0]=wave2stripped[start]
extract/INDX parent,roipnts,FindKruskal(start,parent)==FindKruskal(x,parent)
Make/O/D/N=(numpnts(roipnts)) roi0,roi1,roi2
roi0=wave0stripped[roipnts[x]]
roi1=wave1stripped[roipnts[x]]
roi2=wave2stripped[roipnts[x]]
if(waveexists(wave3stripped))
	Make/O/D/N=(numpnts(roipnts)) roi3
	roi3=wave3stripped[roipnts[x]]
	
endif
end

function/WAVE fetchTree(parent,startPnt)
wave parent
variable startPnt
extract/INDX/FREE parent,res,FindKruskal(startPnt,parent)==FindKruskal(x,parent)
return res
end

function/WAVE quadFromTree(parent,startPnt,dataName)
wave parent
variable startPnt
string dataName
string massSource=dataName+"_0stripped"
string intSource=dataName+"_1stripped"
string timeSource=dataname+"_3stripped"
wave mass=$massSource
wave int=$intSource
Duplicate/O/FREE fetchTree(parent,startPnt) roipnts
variable n=numpnts(roipnts)
make/O/D/N=(n,4)/FREE res
res[][0]=mass[roipnts[x]]
res[][1]=int[roipnts[x]]
res[][2]=defmass(res[x][0])
if(waveexists($timeSource))
	wave time=$timeSource
	res[][3]=time[roipnts[x]]
endif
return res
end

Function formatStds()
ControlInfo/W=GRAPHTTRIBUTOR popup1
string bag=S_value
wave selWave=VariatorSel
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
extract/O/INDX/Free selwave,indexes,selwave==48 && y==0
Make/O/D/N=(114,numpnts(indexes)) Stds_Obs_Mol2add
Stds_Obs_Mol2add=Molbag[x][indexes[y]]
Make/O/D/N=(numpnts(indexes)) Stds_Obs
Stds_Obs=target2theo(ciblebag[indexes[x]][1],ciblebag[indexes[x]][3])
Make/O/D/N=(numpnts(indexes)) Stds_Obs_Charge2Add
Stds_Obs_Charge2Add=chargebag[indexes[x]]
end

Function GrowForest()
wave Stds_Obs
wave/T ResGraphTitle
wave wave1stripped, wave0stripped,wave2stripped
Duplicate/O wave0stripped,Parent_Obs
Parent_Obs=x
wave obs=wave0stripped
variable n=numpnts(obs)
Make/O/D/N=(n,3)/FREE Obs_Vizu
Obs_Vizu[][0]=wave0stripped[x]
Obs_Vizu[][1]=wave1stripped[x]
Obs_Vizu[][2]=wave2stripped[x]
Duplicate/O initiateAdja(Obs) Adja_obs
print "Done Adja"
adja_obs=0
Duplicate/O calcBestof(Obs,Adja_Obs,Stds_Obs) Poids_obs
print "Done Bestof"
Duplicate/O kruskal4mass(obs,Adja_obs,Poids_obs,Parent_Obs,inf) Adja_obs
print "Done Kruskal"
Duplicate/O buildTrace(Adja_obs,Obs_vizu,Parent_Obs), Traceur_obs
print "Done Trace"
Duplicate/O calcDegree(Adja_obs) Degree_Obs
Duplicate/O calcChef(Parent_Obs) Chef_Obs
Duplicate/O CompteTribue(chef_Obs,Adja_Obs,Obs,Parent_Obs) Tribue_Obs
Duplicate/O Sort2D(Tribue_Obs,1,1,-1) Tribue_Obs
Make/O/D/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphSel=0
Make/O/T/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphList=""
ResGraphList=num2str(Tribue_Obs)
end

Function GrowReticles()
wave Stds_Obs
wave/T ResGraphTitle
wave wave1stripped, wave0stripped,wave2stripped
Duplicate/O wave0stripped,Parent_Obs
Parent_Obs=x
wave obs=wave0stripped
variable n=numpnts(obs), m=numpnts(Stds_Obs),k,j
Make/O/D/N=(n,3)/FREE Obs_Vizu
Obs_Vizu[][0]=wave0stripped[x]
Obs_Vizu[][1]=wave1stripped[x]
Obs_Vizu[][2]=wave2stripped[x]
Duplicate/O initiateAdja(Obs) Adja_obs
Duplicate/O/FREE Adja_Obs Adja_temp Poids_temp
print "Done Adja"
Adja_Obs=0
make/O/D/N=1/FREE Stds_Temp
for(j=0;j<m;j+=1)
	Adja_Temp=0
	Stds_Temp[0]=Stds_Obs[j]
	Duplicate/O/FREE calcBestof(Obs,Adja_Temp,Stds_Temp) Poids_temp
	Duplicate/O/FREE kruskal4mass(obs,Adja_Temp,Poids_Temp,Parent_Obs,2) Adja_Temp
	Adja_Obs+=Adja_temp*(j+1)
endfor
Duplicate/O buildTrace(Adja_obs,Obs_vizu,Parent_Obs), Traceur_obs
print "Done Trace"
Duplicate/O calcDegree(Adja_obs) Degree_Obs
Duplicate/O calcChef(Parent_Obs) Chef_Obs
Duplicate/O CompteTribue(chef_Obs,Adja_Obs,Obs,Parent_Obs) Tribue_Obs
Duplicate/O Sort2D(Tribue_Obs,1,1,-1) Tribue_Obs
Make/O/D/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphSel=0
Make/O/T/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphList=""
ResGraphList=num2str(Tribue_Obs)
end

Function GrowGrids()
wave Stds_Obs
wave/T ResGraphTitle
wave wave1stripped, wave0stripped,wave2stripped
Duplicate/O wave0stripped,Parent_Obs
Duplicate/O/FREE Parent_Obs,Parent_temp
Parent_Obs=x
wave obs=wave0stripped
variable n=numpnts(obs), m=numpnts(Stds_Obs),k,j
Make/O/D/N=(n,3)/FREE Obs_Vizu
Obs_Vizu[][0]=wave0stripped[x]
Obs_Vizu[][1]=wave1stripped[x]
Obs_Vizu[][2]=wave2stripped[x]
Duplicate/O initiateAdja(Obs) Adja_obs
Duplicate/O/FREE Adja_Obs Adja_temp Poids_temp
print "Done Adja"
Adja_Obs=0
make/O/D/N=1/FREE Stds_Temp
for(j=0;j<m;j+=1)
	Adja_Temp=0
	Parent_temp=x
	Stds_Temp[0]=Stds_Obs[j]
	Duplicate/O/FREE calcBestof(Obs,Adja_Temp,Stds_Temp) Poids_temp
	Duplicate/O/FREE kruskal4mass(obs,Adja_Temp,Poids_Temp,Parent_temp,2) Adja_Temp
	Adja_Obs+=Adja_temp*(j+1)
	for(k=0;k<n;k+=1)
		if(Parent_obs[k]!=parent_temp[k])
			unionKruskal(Parent_obs[k], parent_temp[k],Parent_obs)
		endif
	endfor
endfor
Duplicate/O buildTrace(Adja_obs,Obs_vizu,Parent_Obs), Traceur_obs
print "Done Trace"
Duplicate/O calcDegree(Adja_obs) Degree_Obs
Duplicate/O calcChef(Parent_Obs) Chef_Obs
Duplicate/O CompteTribue(chef_Obs,Adja_Obs,Obs,Parent_Obs) Tribue_Obs
Duplicate/O Sort2D(Tribue_Obs,1,1,-1) Tribue_Obs
Make/O/D/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphSel=0
Make/O/T/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphList=""
ResGraphList=num2str(Tribue_Obs)
end

function chefs2roi(tribue_obs,mass,int)
wave tribue_obs,mass,int
Make/O/D/N=(dimsize(tribue_obs,0)) roi0,roi1,roi2
roi0=mass[tribue_obs[x][0]]
roi1=int[tribue_obs[x][0]]
roi2=defmass(roi0)
end

function AttribChefs(tribue_obs,mass,int,bag)
wave tribue_obs,mass,int
string bag
Make/O/D/N=(dimsize(tribue_obs,0))/free lesmasses,lesint
lesmasses=mass[tribue_obs[x][0]]
lesint=int[tribue_obs[x][0]]
MolbagAttribList(lesmasses,lesint,bag)
end

function/wave WesslauFitTribues(parent, tribuecrop,mass,int,addresse)//whichone == selected adress of seeds
wave mass,int,parent,tribuecrop
string addresse
nvar charge
variable n=(dimsize(tribuecrop,0)),k,seedmass, norme
Make/O/D/N=(n)/free lesmasses,lesint
lesmasses=mass[tribuecrop[x][0]]
lesint=int[tribuecrop[x][0]]
MolbagKill("temp")
MolbagAttribList(lesmasses,lesint,"temp")
Make/o/d/n=(14+4,n)/free results
wave Molbag_temp,roi0,roi1,ciblebag_temp
wave/T Listbag_temp
Duplicate/o/free canon2weird(Molbag_temp,"Decomposeur") destroyme
results[14,17][]=destroyme[x-14][y]
destroyme[1][]=mod(destroyme[1][y],1)
Duplicate/o/free weird2canon(destroyme,"Decomposeur") destroymetoo
string labelseed, winame
variable intenstotal=sum(int)
NewPath/C/O WesslauPictFolder , addresse
for(k=0;k<n;k+=1)
	extract/o/free destroymetoo,transmol,y==k
	seedmass=mol2massmono(transmol,charge)
	labelseed="\f01"+mol2str(transmol,charge)+" (CH\B2\M)\Bn\M\r\f00"+num2str(ciblebag_temp[k][3])+" ppm error"
	Arbre2ROI(Parent,tribuecrop[k][0])
	Make/o/D/N=5 W_coef
	W_coef[0] = {1e9,2e-8,seedmass,1.03,50}//(tribuecrop[k][2])^1.5}
	norme=norm(roi1)
	roi1/=norme
	FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  roi1 /X=roi0 /D
	FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  roi1 /X=roi0 /D
	FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  roi1 /X=roi0 /D
	results[0,13][k]=condirez()[x]
	results[5][k]=intenstotal
	results[7][k]=sum(roi1)
	labelseed+="\r Residue: "+num2str(results[10][k])
	winame="Fit result "+num2str(k)
	display/k=1/N=winame roi1 vs roi0 as winame
	ModifyGraph mode=8,marker=19,rgb=(34816,34816,34816),useMrkStrokeRGB=1
	SetAxis bottom wavemin(mass),wavemax(mass)
	SetAxis left 0,1.05*wavemax(roi1)
	AppendToGraph fit_roi1
	ModifyGraph lsize(fit_roi1)=2,rgb(fit_roi1)=(0,26112,39168)
	TextBox/C/N=text0/A=RT/X=0.00/Y=0.00 labelseed
	ModifyGraph fStyle=1,axThick=1.5
	Label left "Intensity (a.u.)"
	Label bottom "\\f03m/z\\f01 (u)"
	SavePICT/P=WesslauPictFolder/E=-5/B=288/O as "Fit result "+num2str(k)+".png"
	dowindow/K winame
	//doupdate
endfor
return results
end

function/WAVE WesslauFitsingle(dataY,dataX,bagname)
wave dataX,dataY
string bagname
nvar charge
Make/o/d/n=(14+4,1)/free results
string molbagName="Molbag_"+bagname
wave molbag=$molbagname
variable seedmass,norme
Duplicate/o/free canon2weird(molbag,"Decomposeur") destroyme
results[14,17][]=destroyme[x-14][y]
destroyme[1][]=mod(destroyme[1][y],1)
Duplicate/o/free weird2canon(destroyme,"Decomposeur") destroymetoo
Make/o/D/N=5 W_coef
extract/o/free destroymetoo,transmol,y==0
seedmass=mol2massmono(transmol,charge)
W_coef[0] = {1e9,1e-7,seedmass,1.03,50}
norme=norm(dataY)
dataY/=norme
FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  dataY /X=dataX /D
FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  dataY /X=dataX /D
FuncFit/Q/H="10100"/NTHR=0/TBOX=768 Wesslau W_coef  dataY /X=dataX /D
results[0,13][0]=condirez()[x]
results[5][0]=sum(dataY)
results[7][0]=sum(dataY)
display/k=1 dataY vs dataX
AppendToGraph fit_dataY

return results
end

function WesslauInterface(name,addresse)
string name,addresse
string truename=name
name+="_WesslauFits"
wave/T resgraphlist
wave parent_obs,wave1stripped,wave0stripped,resgraphsel
extract/O/INDX/FREE resgraphsel,whichone,resgraphsel==1
Make/O/D/N=(numpnts(whichone),3)/FREE tribuecrop
tribuecrop=str2num(resgraphlist[whichone[x]][y])
Duplicate/O WesslauFitTribues(parent_obs,tribuecrop ,wave0stripped,wave1stripped,addresse) $name
MolbagKill(truename)
MolbagRename("temp",truename)
end

function wesslauvoir()
string lesfits=WaveList("*_Wesslaufits", ";", "" ),ley,lex
string namebag, name, kulz
variable n=itemsinlist(lesfits),k
variable m,j, indx
display/N=WesslauPlot/K=1
display/N=WesslauError/K=1
make/O/D/N=(0,23) WesslauCrop
make/O/T/N=0 Wesslaulabel
for(k=0;k<n;k+=1)
	name=ReplaceString("_Wesslaufits",stringfromlist(k,lesfits),"")
	namebag="ciblebag_"+name
	kulz="molbag_"+name
	ley=stringfromlist(k,lesfits)
	lex=stringfromlist(k,lesfits)
	wave bag=$namebag
	wave wy=$ley
	wave wx=$lex
	wave molbag=$kulz
	appendtograph/W=WesslauPlot wy[3][*] vs wx[1][*]
	ModifyGraph/W=WesslauPlot mode($ley)=3,marker($ley)=19,useMrkStrokeRGB($ley)=0,zColor($ley)={wy[19][*],*,*,BlueRedGreen,0},zColorMax($ley)=NaN
	//ModifyGraph mode($ley)=3,marker($ley)=19,useMrkStrokeRGB($ley)=0,zColor($ley)={wy[13][*],0,-4,BlueRedGreen,0},zColorMax($ley)=NaN
	ErrorBars/W=WesslauPlot $ley XY,wave=(wy[6][*],wy[6][*]),wave=(wy[8][*],wy[8][*])
	ModifyGraph/W=WesslauPlot log=1
	m=dimsize(wy,1)
	if(dimsize(wy,0)==18)
		insertpoints/M=0 18,5, wy
	endif
	for(j=0;j<m;j+=1)
		if(1)//wy[6][j]<1e-8 && wy[6][j]>1e-10 && wy[8][j]<1e-1 && wy[8][j]>1e-4 && abs(bag[j][3])<5 && wy[13][j]>24)
			indx=dimsize(WesslauCrop,0)
			insertpoints/M=0 indx, 1, WesslauCrop,wesslaulabel
			Wesslaucrop[indx][]=wy[y][j]
			Wesslaucrop[indx][18]=j
			wesslaulabel[indx]=name
		endif
	endfor
	appendtograph/W=WesslauError wy[8][*] vs wx[10][*]
	ModifyGraph/W=WesslauError log=1
endfor
extract/O wesslaucrop, wesslauaux, y==0
display/N=WesslauBulk/K=1 wesslaucrop[][3] vs wesslaucrop[][1]
ModifyGraph/W=WesslauBulk log=1, mode=2, zColor(WesslauCrop)={wesslauaux,*,*,Grays,0}
SetAxis/W=WesslauBulk left 1,1.2
SetAxis/W=WesslauBulk bottom 5.3669045e-10,1.133487e-06
end

function WesslauMAJ()
string lesfits=WaveList("*_Wesslaufits", ";", "" ),ley,lex,name
variable n=itemsinlist(lesfits),k
for(k=0;k<n;k+=1)
	name=ReplaceString("_Wesslaufits",stringfromlist(k,lesfits),"")
	ley=stringfromlist(k,lesfits)
	wave wy=$ley
	ModifyGraph zmrkSize($ley)={wy[13][*],-4,-1,6,1}
endfor
end

function WesslauBackWrite(WesslauLabel, Wesslaucrop)
wave/T Wesslaulabel
wave Wesslaucrop
if(dimsize(wesslaulabel,0)!=dimsize(wesslaucrop,0))
	return 0
endif
variable n=dimsize(wesslaucrop,0),k
string name=""
for(k=0;k<n;k+=1)
	name=Wesslaulabel[k]+"_Wesslaufits"
	wave w=$name
	w[][Wesslaucrop[k][18]]=Wesslaucrop[k][x]
endfor
end

function WesslauTicked(MinNPnts)
variable MinNpnts
variable temps=ticks
string Datalist=TickedDataNameList(),data=""
variable n= itemsinlist(Datalist),k
DoAlert 1,"Have you checked the attribution and graph matrices"
if(V_flag==2)
	return 0
endif
Newpath/C/M="Choose a folder to save pictures"/O/Q WesslauPictFolder
Pathinfo WesslauPictFolder
string rue=S_path
string immeuble
for(k=0;k<n;k+=1)
	data=stringfromlist(k,Datalist)
	immeuble=rue+data+":"
	print immeuble
	activeladata(data)
	GrowForest()
	wave resgraphsel
	wave/T resgraphlist
	resgraphsel=str2num(resgraphlist[x][1])>=MinNPnts && y==0
	WesslauInterface(Data,immeuble)
endfor
temps=(ticks-temps)/60
print num2str(n)+" spectra treated in "+num2str(temps)+" seconds"
end

function WesslauSelected(MinNPnts)
variable MinNpnts
variable temps=ticks
controlinfo/W=advancedmanager TabCon_List_DATA
wave/T datalist=$S_Value
string data=datalist[V_value][1]
DoAlert 1,"Have you checked the attribution and graph matrices"
if(V_flag==2)
	return 0
endif
Newpath/C/M="Choose a folder to save pictures"/O/Q WesslauPictFolder
Pathinfo WesslauPictFolder
string rue=S_path
string immeuble
immeuble=rue+data+":"
print immeuble
activeladata(data)
GrowForest()
wave resgraphsel
wave/T resgraphlist
resgraphsel=str2num(resgraphlist[x][1])>=MinNPnts && y==0
WesslauInterface(Data,immeuble)
temps=(ticks-temps)/60
print "1 spectrum treated in "+num2str(temps)+" seconds"
end

function WesslauActive(Name,MinPnts)
string Name
variable MinPnts
variable temps=ticks
string data=Name
DoAlert 1,"Have you checked the attribution and graph matrices"
if(V_flag==2)
	return 0
endif
Newpath/C/M="Choose a folder to save pictures"/O/Q WesslauPictFolder
Pathinfo WesslauPictFolder
string rue=S_path
string immeuble
immeuble=rue+data+":"
print immeuble
wave resgraphsel
wave/T resgraphlist
resgraphsel=str2num(resgraphlist[x][1])>=MinPnts && y==0
WesslauInterface(Data,immeuble)
temps=(ticks-temps)/60
print "1 spectrum treated in "+num2str(temps)+" seconds"
end

function/S TickedDataNameList()
wave seladvdatalist
wave/T advdatalist
Extract/FREE/T advdatalist,table, y==1 && seladvdatalist[x][5]==48
variable n=dimsize(table,0),k
string liste=""
for(k=0;k<n;k+=1)
	liste+=table[k]+";"
endfor
return liste
end

///////////////////////////////////////////////////////////

Function ButtonProc_17(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=attributeur popup1
			initGRAPHTTRIBUTOR(S_value)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//////////////////////////////////////////////////////

function/WAVE Sort2D(wave2sort,dim,index,sens)
wave wave2sort
variable dim,index,sens
Make/O/D/N=0/FREE res
	switch(dim)
		case 0:
		extract/FREE wave2sort, guide, x==index
			break
		case 1:
		extract/FREE wave2sort, guide, y==index
			break
		default:
			return res
	endswitch
Duplicate/O/FREE guide, waveindex
	if(sens>0)
		makeindex guide, waveindex
	else
		makeindex/R guide, waveindex
	endif
duplicate/O/FREE wave2sort, res
	switch(dim)
		case 0:
		res=wave2sort[x][waveindex[y]]
			break
		case 1:
		res=wave2sort[waveindex[x]][y]
			break
		default:
			return res
	endswitch
return res
end

Function removeHelp()
string wl=winlist("*",";","")
variable n=itemsinlist(wl),k,i,j,test
string lawin,cnl,lacontrol,empty=""
for(k=0;k<n;k+=1)
	lawin=StringFromList(k, wl)
	cnl=ControlNameList(lawin)
	i=itemsinlist(cnl)
	for(j=0;j<i;j+=1)
		lacontrol=StringFromList(j, cnl)
		ControlInfo /W=$lawin $lacontrol
		test=strsearch(s_recreation, " "+lacontrol, 0)
		if(strsearch(s_recreation, "help={", 0) != -1)
			print s_recreation[0,test-1]+" "+lacontrol+" "+"help={\"\"}"+",win="+lawin
			execute s_recreation[0,test-1]+" "+lacontrol+" "+"help={\"\"}"+",win="+lawin
		endif
	endfor
endfor
end

Function wave2molbag(list,name)//transform a text wave into
wave/T list
string name
svar nommol
wave molecule,simu_proba
variable n=dimsize(list,0),k,charge
print n
for(k=0;k<n;k+=1)
	MolbagAddTo(str2mol(list[k]),str2charge(list[k]),mol2str(str2mol(list[k]),str2charge(list[k])),{0,mol2massmono(str2mol(list[k]),str2charge(list[k])),1,0,0},name)
endfor
end

function reshapespecarray(w,wi)
wave w,wi
MatrixOP/O w=w^t
MatrixOP/O wi=wi^t
extract/O w,linw,wi!=0
extract/O wi,linwi,wi!=0
extract/O/INDX linw,steps, linw[x+1]<linw[x]
variable n=numpnts(steps)+1,k
print "nombre de ceo",n
Make/O/D/N=(n) spikes
spikes[0]=steps[0]
for(k=1;k<n-1;k+=1)
	spikes[k]=steps[k]-steps[k-1]
endfor
spikes[n-1]=numpnts(linw)-steps[n-2]-1
variable rows=wavemax(spikes)+1
make/O/D/N=(rows,n) nm,ni
nm=0
ni=0
variable m=numpnts(linw),j=0,l=0,v=-inf
k=0
for(l=0;l<n;l+=1,k+=1)
	for(j=0;j<spikes[l]+1;j+=1,k+=1)
		nm[j][l]=linw[k]
		ni[j][l]=linwi[k]
	endfor
endfor
insertpoints/M=1 0,1,ni,nm
DUplicate/O ni,wi
Duplicate/O nm,w
end

////////////////////////////////////////////////////////////////////////
//                Transforbase
////////////////////////////////////////////////////////////////////////
function/Wave canon2weird(mole,bag)
string bag
wave mole
string smolbag="Molbag_"+bag
wave molbag=$smolbag
Duplicate/O/FREE NonZeroElement(bag) eltz
variable n=numpnts(eltz)
variable m=dimsize(molbag,1)
if(n==m)
	make/O/D/N=(n,n)/FREE w2c
else
	return NaN
endif
w2c=molbag[eltz[x]][y]
variable len=dimsize(mole,1)
make/O/D/N=(n,len) res
res=mole[eltz[x]][y]
MatrixOP/O/FREE res=inv(w2c) x res
return res
end

function/wave weird2canon(mole,bag)
string bag
wave mole
string smolbag="Molbag_"+bag
wave molbag=$smolbag
Duplicate/O/FREE NonZeroElement(bag) eltz
variable n=numpnts(eltz)
variable m=dimsize(molbag,1)
if(1)//n==m)
	make/O/D/N=(n,m)/FREE w2c
else
	return NaN
endif
w2c=molbag[eltz[x]][y]
MatrixOP/O/FREE canon = w2c x mole
Make/O/D/N=(114,dimsize(mole,1))/FREE res
variable k
for(k=0;k<n;k+=1)
	res[eltz[k]][]=canon[k][y]
endfor
return res
end


///////////////////////////////////////////////////////////////////////
//  END  ---  Transforbase
///////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
//			Add to molbag Button
/////////////////////////////////////////////////////////////////////
function/S selectedbags()
string sleslist=wavelist("advmolmanlist*",";","")
string slessel=wavelist("seladvmolman*",";","")
variable n=itemsinlist(sleslist),k,m,j
string res=""
string slalist,slasel
for(k=0;k<n;k+=1)
	slalist=stringfromlist(k,sleslist)
	slasel=stringfromlist(k,slessel)
	wave/T lalist=$slalist
	wave lasel=$slasel
	m=dimsize(lalist,0)
	for(j=0;j<m;j+=1)
		if(lasel[j][0]==65)
			res+=lalist[j][1]+";"
		endif
	endfor
endfor
return res
end

function addmoltomultiplebags(molecule,charge,nom,cible,strlist)
wave molecule,cible
string nom, strlist
variable charge
variable n=itemsinlist(strlist),k
if(n==0)
	strlist="RENAME ME"
	n=1
endif
string lebag
for(k=0;k<n;k+=1)
	lebag=stringfromlist(k,strlist)
	MolbagAddTo(molecule,charge,nom,cible,lebag)
endfor
end

function addcurrentmoltoselected()
wave molecule
nvar charge,lastincox,lastincoy
svar nommol
make/O/D/N=5/FREE cible={0,lastincox,lastincoy,1e6*(mol2massmono(molecule,charge)-lastincox)/mol2massmono(molecule,charge),0}
addmoltomultiplebags(molecule,charge,nommol,cible,selectedbags())
end

Function Button_add_to_selected(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			addcurrentmoltoselected()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

/////////////////////////////////////////////////////////////////////
//		END --- Add to molbag Button
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
//	BEGIN ---- Fasttribution
/////////////////////////////////////////////////////////////////
Function ButtonFasttribution(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			string newname="Fasttribution"
			prompt newname, "Create a new list or add to existing"
			Doprompt "Enter a list name for the newly attributed items to land in",newname
			wave roi0,roi1
			if(v_flag==0)
				MolbagAttribList(roi0,roi1,newname)
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function MolbagAttribList(Masslist,Intlist,bag)
wave Masslist,Intlist
string bag
wave incox, incoy
Variable k, n=numpnts(Masslist),t
Nvar nbprop,inco
nbprop=1
//Make/O/D/N=(n) tempsAttribKillMe=0
for(k=0;k<n;k+=1)
	t=ticks
	inco=Masslist[k]
	incox=Masslist[k]
	incoy=Intlist[k]
	analysemass("arg")
	MolbagExtractor("Propositions",{0},bag)
	//tempsAttribKillMe[k]=(ticks-t)/60.15
endfor
end
/////////////////////////////////////////////////////////////////////
//		END --- Fasttribution
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
//	Plot bag on Agregateur and DMvM
/////////////////////////////////////////////////////////////////
function/WAVE InitTraceMolbag(bag)
string bag
string smolbag="Molbag_"+bag
string schargebag="Chargebag_"+bag
string slistbag="Listbag_"+bag
string sciblebag="Ciblebag_"+bag
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable n=numpnts(chargebag)
MAKE/O/D/N=(n,5)/FREE trace
multithread trace[][0]=target2theo(ciblebag[x][1],ciblebag[x][3])//masse théorique
multithread trace[][1]=ciblebag[x][2]//intensité
multithread trace[][2]=defmass(trace[x][0])//defaut de masse
multithread trace[][3]=abs(ciblebag[x][3])//delta ppm
multithread trace[][4]=ciblebag[x][1]-trace[x][0]//pour barre d'erreur
return trace
end

function OperateTraceOnGraph(bag,win,onoff)
string bag,win
variable onoff
string strace="Traces_"+bag
wave trace=$strace
variable ison
ison=((strsearch(TraceNameList(win, ";", 1 ), strace, 0,2))!=(-1))
	switch(onoff)
		case 1:
			//on la met à jour
			Duplicate/O InitTraceMolbag(bag) $strace
			if(ison==0)
				if(StringMatch(win, "agregateur" ))
					appendtograph/W=$win $strace[][1] vs $strace[][0]
					ModifyGraph/W=$win mode($strace)=8,marker($strace)=60,useMrkStrokeRGB($strace)=1
					ModifyGraph/W=$win zmrkSize($strace)={$strace[*][3],*,*,0,10}
					ErrorBars/W=$win $strace X,wave=($strace[*][4],)
				elseif(StringMatch(win, "dmvm" ))
					appendtograph/W=$win $strace[][2] vs $strace[][0]
					ModifyGraph/W=$win mode($strace)=3,marker($strace)=60,useMrkStrokeRGB($strace)=1
					ModifyGraph/W=$win zmrkSize($strace)={$strace[*][3],*,*,0,10}
				endif
			endif
			break
		case 0:
			if(ison==1)
				Removefromgraph/W=$win $strace
			endif
			break
	endswitch
ison= ((strsearch(TraceNameList("agregateur", ";", 1 ), strace, 0,2))!=(-1)) || ((strsearch(TraceNameList("dmvm", ";", 1 ), strace, 0,2))!=(-1))
if(ison==0)
	Killwaves/Z $strace
endif
end

function/S col2win(col)
variable col
string res
if(col==3)
	res="agregateur"
elseif(col==4)
	res="dmvm"
else
	res=""
endif
return res
end

Function MajTraceOnGraph(selwave,listwave)
wave selwave
wave/T listwave
variable n=dimsize(selwave,0),k
for(k=0;k<n;k+=1)
	OperateTraceOnGraph(listwave[k][1],"agregateur",selwave[k][3]==48)
	OperateTraceOnGraph(listwave[k][1],"dmvm",selwave[k][4]==48)
endfor
end

function EtikTrace(win,prefixe,suffixe,mass,error,onoff)
string win,prefixe,suffixe
variable mass,error,onoff
string lestraces=ListMatch(TraceNameList(win, ";", 1 ),prefixe+"*"+suffixe),latrace
variable n=itemsinlist(lestraces),k,m,j
string bag,letag="",letexte
for(k=0;k<n;k+=1)
	latrace=StringFromList(k,lestraces)
	bag=ReplaceString(prefixe, StringFromList(k,lestraces), "" )
	bag=ReplaceString(suffixe, bag, "" )
	string smolbag="Molbag_"+bag
	string schargebag="Chargebag_"+bag
	string slistbag="Listbag_"+bag
	string sciblebag="Ciblebag_"+bag
	wave molbag=$smolbag
	wave chargebag=$schargebag
	wave/T listbag=$slistbag
	wave ciblebag=$sciblebag
	m=numpnts(chargebag)
	for(j=0;j<m;j+=1)
		letag=bag+"_"+num2str(j)
		if(onoff)
			letexte=listbag[j]
			if(mass)
				letexte+="\r"+num2str(target2theo(ciblebag[k][1],ciblebag[k][3]))
			endif
			if(error)
				letexte+="\r"+num2str(ciblebag[k][3])
			endif
			Tag/W=$win/C/N=$letag/A=LB/Y=(5)/X=(5)/L=2 $latrace , j, letexte
		else
			Tag/W=$win/N=$letag/K
		endif
	endfor
endfor
end

/////////////////////////////////////////////////////////////////
//	END Plot bag on Agregateur
//////////////////////////////////////////////////////////////////


Function Buttoncalculate(ba) : ButtonControl
	STRUCT WMButtonAction &ba
nvar inco, treedepth
wave/Z list_molperm, molecule
string lawave
	switch( ba.eventCode )
		case 2: // mouse up
		if(numpnts(list_molperm)>0)
			analysemass("arg")
			setaxis/W=elaborateur bottom, 0.999*wavemin(simu_mass), 1.001*wavemax(simu_mass)
			setaxis/W=elaborateur left, 0.01*wavemax(simu_proba), wavemax(simu_proba)
			lawave=stringfromlist(0,wavelist("molprop_*",";",""))
			molecule=0
			tuelesspec()
			duplicate/O $lawave molecule
			genestringmol()
			//ListBox list0 selRow=0, win=attributeur //commenté le 13 fevrier 2019 pour obsolescence
			InitMolList(0,"Propositions")
			wave seladvmolman0
			seladvmolman0[][0]=64
			controlinfo/W=dmvm popup0
			generatetree(S_value,treedepth,inco)
		else
			razprop()
		endif
			break
	endswitch

	return 0
End

Window attributeur() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(28,533,448,720)
	ListBox list1,pos={168,81},size={100,106},labelBack=(65535,65535,65535),fSize=10
	ListBox list1,frame=2,listWave=root:list_molperm,mode= 1,selRow= 1
	ListBox list1,special= {0,15,0}
	SetVariable setvar0,pos={3,1},size={336,24},proc=SetVarProcInco,title="\\Z12Current m/z"
	SetVariable setvar0,fSize=16,format="%15.13f",fStyle=1
	SetVariable setvar0,valueBackColor=(51456,44032,58880)
	SetVariable setvar0,limits={-inf,inf,0},value= inco,live= 1
	Button button3,pos={4,112},size={46,20},proc=ButtonRemoveFromPerm,title="Remove"
	Button button3,fColor=(65280,43520,0)
	SetVariable setvar1,pos={241,49},size={174,16},bodyWidth=86,title="Computed results:"
	SetVariable setvar1,limits={0,inf,1},value= nbprop
	ListBox list2,pos={268,81},size={150,106},fSize=10,frame=2,listWave=root:minmax
	ListBox list2,selWave=root:driveminmax,mode= 5,special= {0,15,0}
	Button button5,pos={162,28},size={78,40},proc=Buttoncalculate,title="Calculate"
	Button button5,fColor=(65280,43520,0)
	TitleBox title1,pos={193,67},size={35,13},title="Factors",frame=0
	TitleBox title2,pos={288,67},size={90,13},title="Manual Constraints",frame=0
	Button button0,pos={2,27},size={80,40},proc=ButtonAttribute,title="Attribute"
	Button button0,fColor=(65280,0,0)
	CheckBox checkFindInData,pos={342,5},size={73,14},proc=CheckFindInData_proc,title="Find in data"
	CheckBox checkFindInData,variable= findindata
	Button button6,pos={51,112},size={36,20},proc=ButtonemptyPerm,title="Empty"
	Button button6,fColor=(65280,0,0)
	PopupMenu popup1,pos={241,27},size={177,21},bodyWidth=118,proc=PopMenuProc_2,title="Use this list:"
	PopupMenu popup1,mode=13,popvalue="Classical Organics",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
	Button button1,pos={82,27},size={80,40},proc=ButtonProc_17,title="Graphttribution"
	Button button1,fColor=(52224,52224,52224)
	string molsel="selmol"+num2str(9999)
	string mollist="mollist"+num2str(9999)
	string latitle="moltitles"+num2str(9999)
	ListBox mollist,pos={1,200},size={400,200},labelBack=(47872,47872,47872),proc=PilotMolList//,fSize=10
	ListBox mollist,frame=0,listWave=$mollist,selWave=$molsel
	ListBox mollist,titleWave=$latitle,mode= 9,special= {0,20,0}
	ListBox mollist,widths={0,75,0,30,40,40,0,20,0,0,15},userColumnResize= 1
EndMacro

Function/wave OneNeighbour(data)
wave data
Duplicate/O/FREE data, parent,index
MakeIndex data, index
variable n=numpnts(data),k
for(k=0;k<n;k+=1)
	if((data[index[k+1]]-data[index[k]])<(data[index[k]]-data[index[k-1]]))
		parent[index[k]]=data[index[k+1]]
	else
		parent[index[k]]=data[index[k-1]]
	endif
endfor
//parent=findHK(parent,x)
return parent
end

					  
														  
								 
		   
				 
							  
							  
			 
			 
								   
	  
   


/////////////////////
// Bayesian cut for Kruskal
//////////////////////
function bayesianCut(listepoids)
wave listepoids
variable n=numpnts(listepoids)
variable k,pms,var,avg, stop, c=4,varini
Make/O/D/N=(n) prob=1, altprob=1
varini=variance(listepoids)
//listepoids/=sqrt(varini)
for(k=2;k<n;k+=1)
	var=variance(listepoids,0,k)
	avg=mean(listepoids,0,k)
	//pms=2-2*(0.5*(1+erf((listepoids[k]-avg)/sqrt(2*var))))
	pms=exp(-(listepoids[k]-avg)^2/(2*var))/sqrt(2*pi*var)//gaussian
	//pms=1/(pi*var*(1+((listepoids[k]-avg)/var)^2))//lorentzian
	//pms=sqrt(c/(2*pi))*((exp(-c/(2*(listepoids[k]-avg))))/((listepoids[k]-avg)^(3/2)))//levy distribution
	if(var!=0)
		altprob[k]=pms
		prob[k]=pms*prob[k-1]/(1+2*pms*prob[k-1]-pms-prob[k-1])
		if(prob[k]<0.5)
			stop=k
			k=inf
		endif
	endif
endfor
//listepoids*=sqrt(varini)
return stop
end
/////////////////////
// Bayesian cut for Kruskal
//////////////////////

Window DMVMcnrsCOLORS() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(359.25,73.25,1112.25,671.75) Traceur_obs[*][2] vs Traceur_obs[*][0]
	AppendToGraph simu_probadm vs simu_mass
	AppendToGraph roiy vs roix
	AppendToGraph yplode vs xplode
	AppendToGraph cauchyyplode vs cauchyxplode
	AppendToGraph currentseg[*][0] vs currentseg[*][1]
	AppendToGraph wave2stripped vs wave0stripped
	AppendToGraph roi2 vs roi0
	AppendToGraph defmasstree vs tree
	AppendToGraph incoxdm vs incox
	ModifyGraph gbRGB=(0,10280,19275)
	ModifyGraph mode(simu_probadm)=3,mode(wave2stripped)=3,mode(roi2)=3,mode(defmasstree)=3
	ModifyGraph mode(incoxdm)=3
	ModifyGraph marker(simu_probadm)=18,marker(wave2stripped)=19,marker(roi2)=19,marker(defmasstree)=55
	ModifyGraph marker(incoxdm)=42
	ModifyGraph lSize(yplode)=1.2,lSize(cauchyyplode)=1.2,lSize(currentseg)=2,lSize(roi2)=1.2
	ModifyGraph rgb(Traceur_obs)=(11565,39835,46260),rgb(simu_probadm)=(0,52224,52224)
	ModifyGraph rgb(yplode)=(0,0,65280),rgb(currentseg)=(0,52224,52224),rgb(wave2stripped)=(51400,59110,59110)
	ModifyGraph rgb(roi2)=(36864,14592,58880),rgb(defmasstree)=(34816,34816,34816),rgb(incoxdm)=(26368,0,52224)
	ModifyGraph msize(wave2stripped)=5,msize(roi2)=4,msize(defmasstree)=7,msize(incoxdm)=5
	ModifyGraph mrkThick(wave2stripped)=1,mrkThick(defmasstree)=1
	ModifyGraph opaque(wave2stripped)=1,opaque(incoxdm)=1
	ModifyGraph hideTrace(yplode)=1,hideTrace(cauchyyplode)=1
	ModifyGraph useMrkStrokeRGB(simu_probadm)=1,useMrkStrokeRGB(wave2stripped)=1,useMrkStrokeRGB(roi2)=1
	ModifyGraph useMrkStrokeRGB(defmasstree)=1
	ModifyGraph mrkStrokeRGB(simu_probadm)=(65535,65535,65535),mrkStrokeRGB(wave2stripped)=(0,20560,41120)
	ModifyGraph offset(defmasstree)={-103.374083014535,-0.0879930299186752}
	ModifyGraph zmrkSize(wave2stripped)={logwave1,*,*,1,5},zmrkSize(defmasstree)={treeerror,*,10,6,0}
	ModifyGraph zColor(defmasstree)={treeColor,*,*,Web216,1}
	ModifyGraph zero(left)=1
	ModifyGraph mirror=1
	ModifyGraph fStyle=1
	ModifyGraph axOffset(left)=-2.14286,axOffset(bottom)=-0.388889
	Label left "Mass defect"
	Label bottom "m/z"
	SetAxis left 0.199349310852625,0.414867276670669
	SetAxis bottom 270.330585632119,438.650002513503
	ShowTools/A
	ControlBar 40
	ValDisplay valdisp0,pos={6.00,1.00},size={89.00,17.00},title="ROI size"
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"numpnts(roi1)"
	ValDisplay valdisp1,pos={106.00,1.00},size={99.00,17.00},title="IntCum %"
	ValDisplay valdisp1,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp1,value= #"round(100*(sum(roi1)/sum(wave1stripped)))"
	ValDisplay valdisp2,pos={326.00,1.00},size={74.00,17.00},title="InPnts"
	ValDisplay valdisp2,limits={0,0,0},barmisc={0,1000},value= #"twosig"
	SetVariable setvar0,pos={205.00,1.00},size={120.00,18.00},proc=Proc_setintenscumthres,title="IntCumThres"
	SetVariable setvar0,limits={0,1,0.1},value= intenscumthres
	SetVariable setvar1,pos={6.00,18.00},size={113.00,18.00},proc=SetVarProc_6,title="Depth of tree"
	SetVariable setvar1,limits={0,inf,1},value= _NUM:0
	ValDisplay valdisp3,pos={124.00,18.00},size={142.00,17.00},title="Segment length :"
	ValDisplay valdisp3,limits={0,0,0},barmisc={0,1000}
	ValDisplay valdisp3,value= #"sqrt((currentseg[0][1]-currentseg[1][1])^2+(currentseg[0][0]-currentseg[1][0])^2)"
	PopupMenu popup0,pos={269.00,18.00},size={110.00,19.00},proc=PopMenuProc_3,title="Use this list:"
	PopupMenu popup0,mode=55,popvalue="exva",value= #"replacestring(\"Molbag_\",wavelist(\"Molbag_*\",\";\",\"\"),\"\")"
EndMacro


function rewritecharge(w)
wave/T w
variable n=numpnts(w),k
string olds,news
for(k=0;k<n;k+=1)
	olds=w[k]
	news=replacestring("\M\S1+\M",olds,"\M\S+\M")
	w[k]=news
endfor
end

// Function ADDED CW for modify the RatioPlot + save the resulting image
function RatioColorPerAtoms(molbag,xmin,xmax,ymin,ymax,atom,minatoms,maxatoms,TextDisplay)
string molbag,atom,TextDisplay
variable xmin,xmax,ymin,ymax,minatoms,maxatoms
variable quant=0
string win="AdvMolManOper0"
string sRatiobag=Replacestring(";",ListMatch(TraceNameList(win, ";", 1 ),"*"+molbag+"*"),"")
string sCiblebag="Ciblebag_"+molbag
string sListRatio="ListRatioOper"+num2str(quant)
wave Ratiobag=$sRatiobag
wave Ciblebag=$sCiblebag
wave/T ListRatio=$sListRatio
if(stringmatch(TextDisplay,"Series")==0)
	ListRatio[0][1]=atom+"/C"
	ListRatio[1][1]="Mass"
	InitRatioLists(quant)
endif


extract/FREE ciblebag,transatoms,y==4
make/O/N=(maxatoms-minatoms+1) AtomsNumScale
make/O/T/N=(maxatoms-minatoms+1) AtomsTextScale
AtomsTextScale=num2str(x+1)
AtomsNumScale=(x+1)
ModifyGraph/W=$win axisEnab(bottom)={0,0.9}
SetAxis/W=$win bottom xmin,xmax
SetAxis/W=$win left ymin,ymax
ModifyGraph/W=$win zColor($sRatiobag)={Ciblebag[*][4],minatoms,maxatoms,Red,0}
ColorScale/C/N=text1/F=0/A=MC/W=$win trace=$sRatiobag, axisRange={minatoms,maxatoms}
ColorScale/C/N=text1/W=$win/X=49/Y=0 "Number of "+atom+" atoms"
ColorScale/C/N=text1/W=$win userTicks={AtomsNumScale,AtomsTextScale}
if(stringmatch(TextDisplay,"")==0)
	TextBox/C/N=text0/F=0/A=MC/W=$win/X=47/Y=45 TextDisplay
endif
SavePICT/E=-5/TRAN=1/B=288/M/W=(0,0,14.99,9.984)/Win=$win
end


function RatioColorPerAtomsSeries(molbag,xmin,xmax,ymin,ymax,atom,minatoms,maxatoms,compareatom)
string molbag, atom,compareatom
variable xmin,xmax,minatoms,maxatoms,ymin,ymax
variable k,quant=0
string win="AdvMolManOper0"
string sCiblebag="Ciblebag_"+molbag
wave Ciblebag=$sCiblebag
string sListRatio="ListRatioOper"+num2str(quant)
wave/T ListRatio=$sListRatio
extract/FREE Ciblebag,transatom,y==4
for(k=wavemin(transatom);k<=wavemax(transatom);k+=1)
	if(stringmatch(compareatom,"H"))
		ListRatio[0][1]="H/C*("+atom+"=="+num2str(k)+")"
		ListRatio[1][1]="1+C+N/2-H/2"
	elseif(stringmatch(compareatom,"C"))
		ListRatio[0][1]=atom+"/C*("+atom+"=="+num2str(k)+")"
		ListRatio[1][1]="Mass"
	endif
	InitRatioLists(quant)
	RatioColorPerAtoms(molbag,xmin,xmax,ymin,ymax,atom,minatoms,maxatoms,"Series")
endfor
end


function DBmolbagimportinfo()
print "Valeurs de DBextract"
print "0 : No class"
print "1 : Métabolisme des Acides aminés"
print "2 : Polycétides et peptides particuliers"
print "3 : Métabolites secondaires"
print "4 : Métabolisme des Sucres"
print "5 : Métabolisme énergétique"
print "6 : Métabolisme des Glycan"
print "7 : Métabolisme des lipides"
print "8 : Lipides : Acides gras"
print "9 : Lipides : Glycérolipides"
print "10 : Lipides : Polycétides"
print "11 : Lipides : Prenols"
print "12 : Lipides: Saccharolipides"
print "13 : Lipides: Sphingolipides"
print "14 : Lipides: Sterol lipides"
print "15 : Métabolisme des cofacteurs et vitamines"
print "16 : Métabolisme des nucléotides"
print "17 : Peptides"
print "18 : Di-peptides"
print "19 : Tri-peptides"
print "20 : Tétra-peptides"
print "21 : Xenobiotics Biodegradation and Metabolism"
print "22 : Xenobiotics Drugs etc"
print "23 : 63Mix"
end

function DBmolbagimport(DBextract,chargeion)
variable DBextract, chargeion
string selector
// Détermination de la famille à extraire
if (DBextract==0)
selector=""
elseif(DBextract==1)
selector="Amino Acid Metabolism"
elseif(DBextract==2)
selector="Biosynthesis of Polyketides and Nonribosomal Peptides"
elseif(DBextract==3)
selector="Biosynthesis of Secondary Metabolites"
elseif(DBextract==4)
selector="Carbohydrate Metabolism"
elseif(DBextract==5)
selector="Energy Metabolism"
elseif(DBextract==6)
selector="Glycan Biosynthesis and Metabolism"
elseif(DBextract==7)
selector="Lipid Metabolism"
elseif(DBextract==8)
selector="Lipids: Fatty Acyls"
elseif(DBextract==9)
selector="Lipids: Glycerolipids"
elseif(DBextract==10)
selector="Lipids: Polyketides"
elseif(DBextract==11)
selector="Lipids: Prenols"
elseif(DBextract==12)
selector="Lipids: Saccharolipids"
elseif(DBextract==13)
selector="Lipids: Sphingolipids"
elseif(DBextract==14)
selector="Lipids: Sterol lipids"
elseif(DBextract==15)
selector="Metabolism of Cofactors and Vitamins"
elseif(DBextract==16)
selector="Nucleotide Metabolism"
elseif(DBextract==17)
selector="Peptide"
elseif(DBextract==18)
selector="Peptide(di-)"
elseif(DBextract==19)
selector="Peptide(tri-)"
elseif(DBextract==20)
selector="Peptide(tetra-)"
elseif(DBextract==21)
selector="Xenobiotics Biodegradation and Metabolism"
elseif(DBextract==22)
selector="Xenobiotics Drugs etc"
elseif(DBextract==23)
selector="63Mix"
else
	print "Please put a number between 0 and 23"
	Abort
endif

wave/T/Z CompoundDB_formula, CompoundDB_name, CompoundDB_family//, CompoundDB_iso, CompoundDB_nbiso, CompoundDB_mass 
// Extraction de la famille de composé demandée
if(waveexists(CompoundDB_formula) && waveexists(CompoundDB_name) && waveexists(CompoundDB_family))
	make/O/D/FREE/N=(dimsize(CompoundDB_formula,0)) transselector
	transselector=stringmatch(CompoundDB_family[x],selector)
	// Index des molécules à conserver
	Extract/FREE/INDX transselector, transindx, transselector==1
	// Reconstruction des listes nécessaires à partir des index
	make/O/T/FREE/N=(dimsize(transindx,0)) transformula, transname
	transformula=CompoundDB_formula[transindx[x]]
	transname=CompoundDB_name[transindx[x]]
else
	print "The data base is not loaded or data base have wrong name"
endif

wave molecule
wave/T list_formules
nvar charge
string chargemol, nom
if(chargeion!=0)
	chargemol="H"+num2str(abs(chargeion))
endif
variable k,l=dimsize(transformula,0)
// On nettoye ce qui traine...
videagreg()
// ...et on remplit de nouveau!
for(k=0;k<l;k+=1)
	molecule=str2mol(transformula[k])[x]
	if(chargeion<0)
		molecule-=str2mol(chargemol)[x]
	elseif(chargeion>0)
		molecule+=str2mol(chargemol)[x]
	endif
	charge=chargeion
	tuelesspec()
	genestringmol()
	//addmoltomultiplebags(molecule,charge,nom,cible,"Current")
	addthissimu()
	majliste()
	majlegende()
	list_formules[k]=list_formules[k]+" - "+transname[k]
endfor
nom=UniqueName("DBMolImport",4,0)
fullsavebag(nom)
//videagreg()
end

function isinside(box,vec)
wave box,vec
variable k, n=dimsize(box,0),interior=0
wave w_findlevels
if(n!=0)
	Extract/FREE box,boxy, y==1
	Extract/FREE box,boxx, y==0
	findlevels/Q boxy vec[1]
	duplicate/O/FREE w_findlevels paroix
	paroix=boxx[w_findlevels[x]]-vec[0]<0
	findlevels/Q boxx vec[0]
	duplicate/O/FREE w_findlevels paroiy
	paroiy=boxy[w_findlevels[x]]-vec[1]<0
	//print mod(sum(paroix),2),mod(sum(paroiy),2)
	if(mod(sum(paroix),2)!=0 && mod(sum(paroiy),2)!=0)
		interior+=1
	else
endif
else
endif
return interior
end

function/WAVE side(box,obs)
wave box,obs
variable n=dimsize(obs,0),dim=dimsize(obs,1),k
make/O/D/N=(dim) vec
make/O/D/N=(n) res
for(k=0;k<n;k+=1)
	vec=obs[k][x]
	res[k]=isinside(box,vec)
endfor
return res
end

function/WAVE SpectraFromIonMapPolygons(lin,masse,box)
wave lin,masse,box
variable n=dimsize(lin,0)
Make/O/D/N=(n,2)/FREE obs //masse puis temps
obs=lin[x][y]
Duplicate/O side(box,obs) inout
extract/FREE lin, lesmasses, inout[x] && y==0
extract/FREE lin, lesint, inout[x] && y==3
Duplicate/FREE masse,intensite
intensite=0
variable m=numpnts(lesmasses),k
for(k=0;k<m;k+=1)
	intensite[lesmasses[k]]+=lesint[k]
endfor
return intensite
end



// Expecta-Maxima for chromato

function MainEMchroma(name)
string name
variable t=ticks
string shkwave=name+"_6chroma"
string sleanwave=name+"_12chroma"
wave hkresults=$shkwave
wave leanmap=$sleanwave
Extract hkresults, lesID, y==8
variable n=numpnts(lesID),k
Make/O/FREE/N=(dimsize(leanmap,0)) leanID,lesindex
Multithread leanID=leanmap[x][2]
MakeIndex/O/FREE leanID,lesindex
variable count=0,decount=0
wave EM_Classes, EM_Weights, EM_Cov
for(k=0;k<n;k+=1)
	for(count=count;leanmap[lesindex[count]][2]!=lesID[k];count+=1)
	endfor
	decount=hkresults[k][0]
	Make/O/I/N=(decount,2) Obs2kill
	Obs2kill=leanmap[lesindex[count+x]][y]
	Make/O/I/N=(decount,1) Inten2kill
	Inten2kill=leanmap[lesindex[count+x]][3]
	formaObs(obs2kill,inten2kill)
	if(dimsize(EM_Classes,2)>1)
		//print hkresults[k][0],hkresults[k][1],hkresults[k][2],hkresults[k][3],hkresults[k][4],hkresults[k][5],hkresults[k][6],hkresults[k][7],hkresults[k][8]
		cycle2D(Obs2kill,Inten2kill,EM_Classes,EM_Weights,EM_Cov)
		ModIslandList(name,colorclasses,k,lesID[k],Obs2kill,Inten2kill,lesindex,count,decount)
		//k=n
	endif
	count+=decount
endfor
print "Expecta-Maxima processed in ",(ticks-t)/60,"seconds"
//listelesdonnees()
end

function formaObs(obs2kill,inten2kill)
wave obs2kill,inten2kill
matrixOP/O/FREE minz=mincols(Obs2kill)
matrixOP/O/FREE maxz=maxcols(Obs2kill)
Obs2kill[][1]*=(maxz[0]-minz[0])/(maxz[1]-minz[1])
variable nbin=maxz[0]-minz[0]
Make/O/I/FREE/N=(nbin+1) lesy=0
variable k,n=dimsize(obs2kill,0)
for(k=0;k<n;k+=1)
	lesy[obs2kill[k][0]-minz[0]]+=inten2kill[k]
endfor
extract/O/FREE/INDX lesy,res, lesy[x-1]<lesy[x] && lesy[x]>lesy[x+1]
Setclasses2D(numpnts(res),Obs2kill,inten2kill)
wave EM_classes
EM_classes[0][0][]=res[z]+minz[0]
EM_classes[0][1][]=((maxz[1]+minz[1])/2)*(maxz[0]-minz[0])/(maxz[1]-minz[1])
//inten2kill=1
end

function Setclasses2D(nb,obs,inten)
variable nb
wave obs,inten
variable n=dimsize(obs,0),dim=dimsize(obs,1)
Make/O/D/N=(1,dim,nb) EM_Classes
Make/O/D/N=(n,1,nb) EM_Weights=1
MatrixOP/O/FREE tcov= (obs)^t x (obs)
Make/O/D/N=(dim,dim,nb) EM_Cov
EM_Cov=(tcov[x][y])/sum(inten)
EM_Cov=(x==y)*1+2*x+0.001
end

function cycle2D(obs,inten,EM_Classes,EM_Weights,EM_Cov)
wave obs,inten,EM_Classes,EM_Weights,EM_Cov
variable k,n=400
variable dd=inf
make/O/D/N=(n) likelihood=inf
for(k=0;k<n && dd!=0;k+=1)
	dd=likelihood[max(0,k-1)]
	EM2D(obs,inten,EM_Classes,EM_Weights,EM_Cov)
	likelihood[k]=norm(em_weights)^2/dimsize(obs,0)
	dd-=likelihood[k]
	//doupdate
endfor
colorize2D(EM_weights)
end

function EM2D(pos,inten,EM_Classes,EM_Weights,EM_Cov)
wave pos,inten,EM_Classes,EM_Weights,EM_Cov
variable n=dimsize(pos,0),dim=2, nb=dimsize(EM_Classes,2)
Make/O/D/FREE/N=(n,dim,nb) vecred,weiobs
Make/O/D/FREE/N=(dim,dim,nb) Pass=1,Diag=0,invDiag=0
Make/O/D/FREE/N=(1,1,nb) canon
canon=sqrt(EM_Cov[0][0][r]^2 - 2*EM_Cov[0][0][r]*EM_Cov[1][1][r] +4*EM_Cov[0][1][r]^2 +EM_Cov[1][1][r]^2)
diag[0][0][]=0.5*(-canon[0][0][r]+EM_Cov[0][0][r]+EM_Cov[1][1][r])
diag[1][1][]=0.5*(canon[0][0][r]+EM_Cov[0][0][r]+EM_Cov[1][1][r])
pass[0][0][]=-(-EM_Cov[0][0][r]+EM_Cov[1][1][r]+canon[0][0][r])/(2*EM_Cov[0][1][r])
pass[0][1][]=-(-EM_Cov[0][0][r]+EM_Cov[1][1][r]-canon[0][0][r])/(2*EM_Cov[0][1][r])
MatrixOP/O/FREE transpass=normalizecols(pass)
pass=transpass
diag=sign(diag)*sqrt(abs(diag))
invDiag[0][0][]=1/Diag[0][0][r]
invDiag[1][1][]=1/Diag[1][1][r]
MatrixOP/O/FREE invA= pass x invDiag x inv(pass)
MatrixOP/O/FREE normCov=sqrt(sumcols(sumrows(powR(EM_Cov,2))))
Multithread vecred=pos[x][y]-EM_Classes[0][y][r]
MatrixOP/O/FREE sqmod=sumrows(powR(vecred x invA,2))
Multithread EM_Weights=normCov[0][0][r]^(-0.5)*exp(-0.5*sqmod[x][0][r])/(2*pi)
MatrixOP/O/FREE transEM_Weights=powR(abs(EM_Weights),2)
MatrixOP/O/FREE tirlibibi=sumBeams(transEM_Weights)
MatrixOP/O/FREE transEM_Weights=(EM_Weights)/sqrt(tirlibibi)
Multithread EM_Weights=inten[x]*transEM_Weights[x][0][r]
Multithread weiobs=pos[x][y]*EM_Weights[x][0][r]
MatrixOP/O/FREE lespoids=sumcols(EM_Weights)
MatrixOP/O EM_Classes=sumcols(weiobs)
Multithread EM_Classes/=lespoids[0][0][r]
Multithread vecred=(pos[x][y]-EM_Classes[0][y][r])*sqrt(EM_Weights[x][0][r])
MatrixOP/O EM_Cov=vecred^t x vecred
Multithread EM_Cov/=lespoids[0][0][r]
end

Function colorize2D(EM_weights)
wave EM_weights
Duplicate/O/FREE EM_weights tcolor
variable n=dimsize(EM_weights,0),k
variable nb=dimsize(EM_weights,2),xeu
Make/O/D/FREE/N=(nb) transpoids
for(k=0;k<n;k+=1)
	transpoids=EM_weights[k][0][x]
	xeu=wavemax(transpoids)
	tcolor[k][0][]=(EM_weights==xeu)*r
endfor
MatrixOP/O colorClasses=sumbeams(tcolor)
end

function ModIslandList(name,colorclasses,k,IslandID,Obs2kill,Inten2kill,lesindex,count,decount)
string name
wave colorclasses, obs2kill,inten2kill,lesindex
variable k, islandID,count,decount
//On redresse la dimension temporelle ajustée par EM précédement
matrixOP/O/FREE minz=mincols(Obs2kill)
matrixOP/O/FREE maxz=maxcols(Obs2kill)
Obs2kill[][1]/=(maxz[0]-minz[0])/(maxz[1]-minz[1])
//On appelle ce qu'il faut
string smass2indx=name+"_4chroma"
wave mass2indx=$smass2indx //mass <-> mass index
// Décalage temporel due au offset
string sTime=name+"_3chroma" 
string sMap=name+"_5chroma" 
wave chromatime=$sTime
wave chromamap=$sMap 
variable timeoffset=dimoffset(chromamap,1)
findlevel/Q chromatime, timeoffset
variable timeoffsetindx=ceil(V_levelX)
//******************//
//On modifie 6chroma
string shklist=name+"_6chroma"
wave hklist=$shklist
//on modifie la ligne déjà existante...
Extract/O/INDX/FREE colorClasses,transx,colorClasses==0
Duplicate/O/FREE transx,transval
hklist[k][0]=numpnts(transx)
transval=inten2kill[transx[x]]
Extract/O/FREE/INDX inten2kill,maxintindx,inten2kill==wavemax(inten2kill)
hklist[k][1]=sum(transval)
transval=Obs2kill[transx[x]][0]
hklist[k][2]=mass2indx[wavemin(transval)]
hklist[k][3]=mass2indx[wavemax(transval)]
hklist[k][6]=(hklist[k][2]+hklist[k][3])/2
transval=Obs2kill[transx[x]][1]
hklist[k][4]=chromatime[wavemin(transval)+timeoffsetindx]
hklist[k][5]=chromatime[wavemax(transval)+timeoffsetindx]
hklist[k][7]=chromatime[Obs2kill[maxintindx[0]][1]+timeoffsetindx]
hklist[k][8]=hklist[k][8]
hklist[k][9]=defmass(hklist[k][6])
//... et on ajoute les autres à la fin
variable n=wavemax(colorclasses)
variable m,dimhk=dimsize(hklist,0) // !! Ne pas changer la valeur de dimhk !!
insertpoints/M=0 dimhk,n,hklist
for(m=0;m<n;m+=1)
	Extract/O/INDX/FREE colorClasses,transx,colorClasses==m+1
	Duplicate/O/FREE transx,transval
	hklist[m+dimhk][0]=numpnts(transx)
	transval=inten2kill[transx[x]]
	Extract/O/FREE/INDX inten2kill,maxintindx,inten2kill==wavemax(inten2kill)
	hklist[m+dimhk][1]=sum(transval)
	transval=Obs2kill[transx[x]][0]
	hklist[m+dimhk][2]=mass2indx[wavemin(transval)]
	hklist[m+dimhk][3]=mass2indx[wavemax(transval)]
	hklist[m+dimhk][6]=(hklist[m+dimhk][2]+hklist[m+dimhk][3])/2
	transval=Obs2kill[transx[x]][1]
	hklist[m+dimhk][4]=chromatime[wavemin(transval)+timeoffsetindx]
	hklist[m+dimhk][5]=chromatime[wavemax(transval)+timeoffsetindx]
	hklist[m+dimhk][7]=chromatime[Obs2kill[maxintindx[0]][1]+timeoffsetindx]
	hklist[m+dimhk][8]=hklist[dimhk-1][8]+m+1
	hklist[m+dimhk][9]=defmass(hklist[m+dimhk][6])
endfor
//******************//
//On modifie 12chroma
string sleanmap=name+"_12chroma"
wave leanmap=$sleanmap
variable l
//print "dimhk",dimhk,count,decount
for(l=0;l<decount;l+=1)
	//if(mod(l,500)==0)
		//print l,leanmap[lesindex[count+l]][2]
		//print leanmap[lesindex[count+l]][2],colorclasses[l]
		//print dimhk,hklist[dimhk-1][8],colorclasses[l]
	//endif
	leanmap[lesindex[count+l]][2]=((leanmap[lesindex[count+l]][2]*(colorclasses[l]==0)) + ((hklist[dimhk-1][8]+colorclasses[l])*(colorclasses[l]!=0)))
endfor
end

// Fonctions to compare 2 molbags stoichiometric formulas
								  
// Testing purpose - Box selection for ion map
Window ChromaBoxSelector() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(302,804,697,1012) C18MeOH_12chroma[*][1] vs C18MeOH_12chroma[*][0]
	AppendToGraph box[*][1] vs box[*][0]
	ModifyGraph mode(C18MeOH_12chroma)=3
	ModifyGraph marker(C18MeOH_12chroma)=19
	ModifyGraph rgb(box)=(26411,1,52428)
	ModifyGraph msize(C18MeOH_12chroma)=1
	ModifyGraph zColor(C18MeOH_12chroma)={inout,0,1,YellowHot}
	ModifyGraph axOffset(bottom)=-0.0769231
EndMacro

function CompareMolBag(sref,scomp)
string sref,scomp
string slistref="Listbag_"+sref
string schargeref="Chargebag_"+sref
string scibleref="Ciblebag_"+sref
string smolref="Molbag_"+sref
string slistcomp="Listbag_"+scomp
string schargecomp="Chargebag_"+scomp
string sciblecomp="Ciblebag_"+scomp
string smolcomp="Molbag_"+scomp
wave/T listref=$slistref
wave/T listcomp=$slistcomp
wave chargeref=$schargeref
wave chargecomp=$schargecomp
wave cibleref=$scibleref
wave ciblecomp=$sciblecomp
wave molref=$smolref
wave molcomp=$smolcomp
end

function CompEnsemble(nom1,nom2)
string nom1, nom2
string smol1="Molbag_"+nom1
string smol2="Molbag_"+nom2
wave mol1=$smol1
wave mol2=$smol2
variable n=dimsize(mol1,1)
variable m=dimsize(mol2,1)
variable k=dimsize(mol1,0)
Make/O/D/N=(n,m,k)/FREE compicomp
Multithread compicomp=(mol1[z][x]-mol2[z][y])^2
MatrixOP/O pouet=sumbeams(compicomp)
MatrixOP/O res2=minCols(pouet)
MatrixOP/O res1=minRows(pouet)
Extract/O/INDX res1,communs,res1==0
Extract/O/INDX res1,dans1,res1!=0
Extract/O/INDX res2,dans2,res2!=0
MolbagKill("Common")
MolbagKill("In1Only")
MolbagKill("In2Only")
MolbagExtractor(nom1,communs,"Common")
MolbagExtractor(nom1,dans1,"In1Only")
MolbagExtractor(nom2,dans2,"In2Only")
end

///////////////////////////
//// A new Hough
///////////////////////////
///Interface
function NewHoughInit()
end

function NewHoughSetcontrols()
end

function NewHoughHook()
end

function NewHoughMaj()
end

function NewHoughCleanUp()
end

///Metier

function initNewHoughWaves()
end

function/WAVE FastLinSloInt(masse)
wave masse
variable n=dimsize(masse,0)
make/o/d/n=((n^2-n)/2)/free defm,slopes,massk,massj,defk,defj,interc
MAKE/O/D/N=((n^2-n)/2,5)/FREE res
multithread res[][2] = (floor(0.5*(1+sqrt(8*x+1))))
multithread res[][3] = x-((floor(0.5*(1+sqrt(8*x+1))))^2-(floor(0.5*(1+sqrt(8*x+1)))))/2
multithread massk=masse[res[x][2]][y]
multithread massj=masse[res[x][3]][y]
multithread defk=defmass(massk)
multithread defj=defmass(massj)
FastOP defm = defk - defj
FastOP slopes = massk - massj
multithread res[][4] = slopes[x]
FastOP slopes = defm / slopes
FastOP interc = defj - slopes * massj
multithread res[][0] = slopes[x]
multithread res[][1] = interc[x]
return res
end

function/WAVE NewHoughTransform(masse,minmax,reso)
wave masse,minmax
variable reso
Duplicate/O FastLinSloInt(masse) houghdata
MAKE/O/D/N=(reso,reso)/FREE HoughMap=0
setscale x,minmax[0],minmax[1], HoughMap
setscale y,minmax[2],minmax[3], HoughMap
variable n=dimsize(houghdata,0),k, xx,yy
for(k=0;k<n;k+=1)
	yy=ScaleToIndex(HoughMap,houghdata[k][1],1)
	xx=ScaleToIndex(HoughMap,houghdata[k][0],0)
	if(xx>=0)
		if(xx<reso)
			if(yy>=0)
				if(yy<reso)
					HoughMap[xx][yy]+=1
				endif
			endif
		endif
	endif
endfor
return HoughMap
end

function NewHoughSieve(nomdebag,prof)
wave prof
string nomdebag
string Smolbag="Molbag_"+nomdebag
string Schargebag="Chargebag_"+nomdebag
string elem
wave molbag=$Smolbag
wave chargebag=$Schargebag
wave molecule, mendev_masses
wave/T elements
variable n=dimsize(molbag,1),k,j
Duplicate/O/FREE prof, cote,stretch
stretch=1
cote=(2*prof[x]+1)
for(j=1;j<n;j+=1)
	stretch[j]=stretch[j-1]*cote[j-1]
endfor
variable taille=stretch[n-1]*(cote[n-1])
nvar charge
Make/O/D/N=(n)/FREE lesmasses
Duplicate/O/FREE molecule transmolecule
for(k=0;k<n;k+=1)
	transmolecule=molbag[x][k]*mendev_masses[x][0]/(max(1,abs(chargebag[x])))
	lesmasses[k]=sum(transmolecule)
endfor
Make/O/D/N=(taille,n)/FREE cribletree
Multithread cribletree=(mod(floor(x/(stretch[y])),cote[y])-prof[y])*lesmasses[y]
MatrixOP/O/FREE tree=sumrows(cribletree)
Duplicate/O/FREE tree, defmasstree
Multithread defmasstree=defmass(tree)
Make/O/D/N=(114,taille) Molbag_Hough=0
for(k=0;k<n;k+=1)
	Molbag_Hough+=(mod(floor(y/(stretch[k])),cote[k])-prof[k])*molbag[x][k]
endfor
Make/O/D/N=(taille,5) Ciblebag_Hough=0
Multithread Ciblebag_Hough[][1]=tree[x]
Multithread Ciblebag_Hough[][2]=exp(-tree[x])
Multithread Ciblebag_Hough[][3]=0
Multithread Ciblebag_Hough[][4]=defmasstree[x]/tree[x]
Make/O/D/N=(taille) Chargebag_Hough=charge
Make/O/D/T/N=(taille) Listbag_Hough=""
for(j=0;j<taille;j+=1)
	for(k=0;k<114;k+=1)//verifie la presence de chaque element
			if (Molbag_Hough[k][j]!=0) // si l' element k est present
				elem=elements[k]
				Listbag_Hough[j]+=elem+"\B"+num2str(Molbag_Hough[k][j])+"\M"//genere la formule brute
			endif
		endfor
		if (charge>0)
			Listbag_Hough[j]+="\S"+num2str(abs(charge))+"+\M"
		elseif (charge<0)
			Listbag_Hough[j]+="\S"+num2str(abs(charge))+"-\M"
		endif
endfor
Extract/O/INDX/FREE Ciblebag_Hough, targetmasses, y==1 && Ciblebag_Hough<=0
targetmasses-=taille
MolbagWithdrawFrom(targetmasses,"Hough")
MolbagResort("Hough",7)
taille=dimsize(Ciblebag_Hough,0)
Make/O/D/N=114/FREE mol1,mol2
Make/O/D/N=(taille)/FREE linetokill=0
variable score,mark=0
for(j=1;j<taille;j+=1)
	mol2=Molbag_Hough[x][j]
	mol1=Molbag_Hough[x][j-1]
	score=MatrixDot(mol1,mol2)/(norm(mol2)*norm(mol1))
	if(abs(score-1)<1e-14)
		linetokill[j]=1
		if(norm(mol2)<norm(mol1))
			Ciblebag_Hough[mark][]=Ciblebag_Hough[j][y]
			Listbag_Hough[mark]=Listbag_Hough[j]
			Chargebag_Hough[mark]=Chargebag_Hough[j]
			Molbag_Hough[][mark]=Molbag_Hough[x][j]
		endif
	else
		mark=j
	endif
endfor
Extract/O/INDX/FREE linetokill, targetmasses, linetokill[x]==1
MolbagWithdrawFrom(targetmasses,"Hough")
end


Function MolbagSlope2fifth(name)
string name
string smolbag="Molbag_"+name
string schargebag="Chargebag_"+name
string slistbag="Listbag_"+name
string sciblebag="Ciblebag_"+name
wave molbag=$smolbag
wave chargebag=$schargebag
wave/T listbag=$slistbag
wave ciblebag=$sciblebag
variable k,n=dimsize(ciblebag,0),masse
for(k=0;k<n;k+=1)
	masse=target2theo(ciblebag[k][1],ciblebag[k][3])//mol2massmono(molecule,chargebag[k])//
	ciblebag[k][4]=defmass(masse)/masse
endfor
end


///////////////////////////
//// END   ---   A new Hough
///////////////////////////


//////////////////
/////	Partitioning
//////////////////

//les sets sont les ensembles
//le set a une taille et un ordre
//les elements sont les éléments dans ces ensembles

//un test binaire est un test sur un éléments qui renvoi 0 ou 1

function initPartit(quant)
variable quant

end

////////////
///// End --- Partitioning
///////////


Window WesslauPlot() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(397.5,105.5,1588.5,1077.5)/K=1  Toxi2_ROIed_ROIed_WesslauFits[3][*] vs Toxi2_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph Toxi3_ROIed_ROIed_WesslauFits[3][*] vs Toxi3_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph Toxi2_crop_ROIed_ROI_WesslauFits[3][*] vs Toxi2_crop_ROIed_ROI_WesslauFits[1][*]
	AppendToGraph TOXI1_ROIed_ROIed_RO_WesslauFits[3][*] vs TOXI1_ROIed_ROIed_RO_WesslauFits[1][*]
	AppendToGraph MZUVirrLEGIT_ROIed_WesslauFits[3][*] vs MZUVirrLEGIT_ROIed_WesslauFits[1][*]
	AppendToGraph LS_ROIed_WesslauFits[3][*] vs LS_ROIed_WesslauFits[1][*]
	AppendToGraph OrgueilLDI_WesslauFits[3][*] vs OrgueilLDI_WesslauFits[1][*]
	AppendToGraph JHU112216_ROIed_WesslauFits[3][*] vs JHU112216_ROIed_WesslauFits[1][*]
	AppendToGraph Ivuna50_500_ROIed_RO_WesslauFits[3][*] vs Ivuna50_500_ROIed_RO_WesslauFits[1][*]
	AppendToGraph PicachuP072A_ROIed_WesslauFits[3][*] vs PicachuP072A_ROIed_WesslauFits[1][*]
	AppendToGraph GRA_ROIed_WesslauFits[3][*] vs GRA_ROIed_WesslauFits[1][*]
	AppendToGraph PicachuP072B_ROIed_WesslauFits[3][*] vs PicachuP072B_ROIed_WesslauFits[1][*]
	AppendToGraph Sapp4_ROIed_ROIed_WesslauFits[3][*] vs Sapp4_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph CBokk50_500_ROIed_RO_WesslauFits[3][*] vs CBokk50_500_ROIed_RO_WesslauFits[1][*]
	AppendToGraph JHU112916_ROIed_WesslauFits[3][*] vs JHU112916_ROIed_WesslauFits[1][*]
	AppendToGraph hortholin_ROIed_ROIe_WesslauFits[3][*] vs hortholin_ROIed_ROIe_WesslauFits[1][*]
	AppendToGraph Vass_ROIed_ROIed_ROI_WesslauFits[3][*] vs Vass_ROIed_ROIed_ROI_WesslauFits[1][*]
	AppendToGraph MET_ROIed_WesslauFits[3][*] vs MET_ROIed_WesslauFits[1][*]
	AppendToGraph psocmeto_uf_cal_WesslauFits[3][*] vs psocmeto_uf_cal_WesslauFits[1][*]
	AppendToGraph SappTEMrebuild_crop_WesslauFits[3][*] vs SappTEMrebuild_crop_WesslauFits[1][*]
	AppendToGraph EET_ROIed_WesslauFits[3][*] vs EET_ROIed_WesslauFits[1][*]
	AppendToGraph thoSflow_ROIed_ROIed_WesslauFits[3][*] vs thoSflow_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph METmetolu2_ROIed_WesslauFits[3][*] vs METmetolu2_ROIed_WesslauFits[1][*]
	AppendToGraph JHU102516_ROIed_WesslauFits[3][*] vs JHU102516_ROIed_WesslauFits[1][*]
	AppendToGraph TL5b_ROIed_ROIed_ROI_WesslauFits[3][*] vs TL5b_ROIed_ROIed_ROI_WesslauFits[1][*]
	AppendToGraph GRAmetolu_ROIed_WesslauFits[3][*] vs GRAmetolu_ROIed_WesslauFits[1][*]
	AppendToGraph QUE_ROIed_WesslauFits[3][*] vs QUE_ROIed_WesslauFits[1][*]
	AppendToGraph Orgueil2018_ROIed_Cr_WesslauFits[3][*] vs Orgueil2018_ROIed_Cr_WesslauFits[1][*]
	AppendToGraph JHU122016_ROIed_WesslauFits[3][*] vs JHU122016_ROIed_WesslauFits[1][*]
	AppendToGraph EETmetolu_ROIed_ROIe_WesslauFits[3][*] vs EETmetolu_ROIed_ROIe_WesslauFits[1][*]
	AppendToGraph QUEmetolu2_ROIed_WesslauFits[3][*] vs QUEmetolu2_ROIed_WesslauFits[1][*]
	AppendToGraph QUEmetolu_ROIed_WesslauFits[3][*] vs QUEmetolu_ROIed_WesslauFits[1][*]
	AppendToGraph Renazzo50_500_ROIed__WesslauFits[3][*] vs Renazzo50_500_ROIed__WesslauFits[1][*]
	AppendToGraph MaiaCONGR34_ROIed_WesslauFits[3][*] vs MaiaCONGR34_ROIed_WesslauFits[1][*]
	AppendToGraph TL11i_ROIed_ROIed_cr_WesslauFits[3][*] vs TL11i_ROIed_ROIed_cr_WesslauFits[1][*]
	AppendToGraph JHU120616_ROIed_ROIe_WesslauFits[3][*] vs JHU120616_ROIed_ROIe_WesslauFits[1][*]
	AppendToGraph METmetolu_ROIed_WesslauFits[3][*] vs METmetolu_ROIed_WesslauFits[1][*]
	AppendToGraph NMP20_ROIed_WesslauFits[3][*] vs NMP20_ROIed_WesslauFits[1][*]
	AppendToGraph Murch_Neg_UF_ROIed_N_WesslauFits[3][*] vs Murch_Neg_UF_ROIed_N_WesslauFits[1][*]
	AppendToGraph Fre211b_ROIed_ROIed_WesslauFits[3][*] vs Fre211b_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph EETmetolu2_ROIed_WesslauFits[3][*] vs EETmetolu2_ROIed_WesslauFits[1][*]
	AppendToGraph MaiaCO28_ROIed_WesslauFits[3][*] vs MaiaCO28_ROIed_WesslauFits[1][*]
	AppendToGraph NMP5_ROIed_WesslauFits[3][*] vs NMP5_ROIed_WesslauFits[1][*]
	AppendToGraph thoDflow_ROIed_ROIed_WesslauFits[3][*] vs thoDflow_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph GRAmetolu2_ROIed_WesslauFits[3][*] vs GRAmetolu2_ROIed_WesslauFits[1][*]
	AppendToGraph MACR5_ROIed_WesslauFits[3][*] vs MACR5_ROIed_WesslauFits[1][*]
	AppendToGraph JunkoIsa2_ROIed_ROIe_WesslauFits[3][*] vs JunkoIsa2_ROIed_ROIe_WesslauFits[1][*]
	AppendToGraph NWApos_UF_ROIed_WesslauFits[3][*] vs NWApos_UF_ROIed_WesslauFits[1][*]
	AppendToGraph MACR10_ROIed_WesslauFits[3][*] vs MACR10_ROIed_WesslauFits[1][*]
	AppendToGraph MaiaCOGR30_ROIed_WesslauFits[3][*] vs MaiaCOGR30_ROIed_WesslauFits[1][*]
	AppendToGraph Murch_clean_ROIed_WesslauFits[3][*] vs Murch_clean_ROIed_WesslauFits[1][*]
	AppendToGraph MACR15_ROIed_ROIed_WesslauFits[3][*] vs MACR15_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph Fre311b_ROIed_WesslauFits[3][*] vs Fre311b_ROIed_WesslauFits[1][*]
	AppendToGraph MACR20_ROIed_ROIed_WesslauFits[3][*] vs MACR20_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph MVK5_ROIed_WesslauFits[3][*] vs MVK5_ROIed_WesslauFits[1][*]
	AppendToGraph Fre1011b_ROIed_WesslauFits[3][*] vs Fre1011b_ROIed_WesslauFits[1][*]
	AppendToGraph MVK0_ROIed_WesslauFits[3][*] vs MVK0_ROIed_WesslauFits[1][*]
	AppendToGraph Isono61WN_ROIed_WesslauFits[3][*] vs Isono61WN_ROIed_WesslauFits[1][*]
	AppendToGraph Nogoya50_500_ROIed_R_WesslauFits[3][*] vs Nogoya50_500_ROIed_R_WesslauFits[1][*]
	AppendToGraph fresneau_ROIed_WesslauFits[3][*] vs fresneau_ROIed_WesslauFits[1][*]
	AppendToGraph SucreHCL_ROIed_recal_WesslauFits[3][*] vs SucreHCL_ROIed_recal_WesslauFits[1][*]
	AppendToGraph Fre315b_ROIed_WesslauFits[3][*] vs Fre315b_ROIed_WesslauFits[1][*]
	AppendToGraph NMP0_ROIed_WesslauFits[3][*] vs NMP0_ROIed_WesslauFits[1][*]
	AppendToGraph Isono64WON_ROIed_WesslauFits[3][*] vs Isono64WON_ROIed_WesslauFits[1][*]
	AppendToGraph Fre211bHYDRO_ROIed_R_WesslauFits[3][*] vs Fre211bHYDRO_ROIed_R_WesslauFits[1][*]
	AppendToGraph MACR0_ROIed_WesslauFits[3][*] vs MACR0_ROIed_WesslauFits[1][*]
	AppendToGraph steack_ROIed_ROIed_WesslauFits[3][*] vs steack_ROIed_ROIed_WesslauFits[1][*]
	AppendToGraph PSOCLDI_ROIed_WesslauFits[3][*] vs PSOCLDI_ROIed_WesslauFits[1][*]
	AppendToGraph MZUVtemLEGIT_ROIed_R_WesslauFits[3][*] vs MZUVtemLEGIT_ROIed_R_WesslauFits[1][*]
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph msize(Toxi2_ROIed_ROIed_WesslauFits)=2,msize(Toxi3_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(Toxi2_crop_ROIed_ROI_WesslauFits)=2,msize(TOXI1_ROIed_ROIed_RO_WesslauFits)=2
	ModifyGraph msize(MZUVirrLEGIT_ROIed_WesslauFits)=2,msize(LS_ROIed_WesslauFits)=2
	ModifyGraph msize(OrgueilLDI_WesslauFits)=2,msize(JHU112216_ROIed_WesslauFits)=2
	ModifyGraph msize(Ivuna50_500_ROIed_RO_WesslauFits)=2,msize(PicachuP072A_ROIed_WesslauFits)=2
	ModifyGraph msize(GRA_ROIed_WesslauFits)=2,msize(PicachuP072B_ROIed_WesslauFits)=2
	ModifyGraph msize(Sapp4_ROIed_ROIed_WesslauFits)=2,msize(CBokk50_500_ROIed_RO_WesslauFits)=2
	ModifyGraph msize(JHU112916_ROIed_WesslauFits)=2,msize(hortholin_ROIed_ROIe_WesslauFits)=2
	ModifyGraph msize(Vass_ROIed_ROIed_ROI_WesslauFits)=2,msize(MET_ROIed_WesslauFits)=2
	ModifyGraph msize(psocmeto_uf_cal_WesslauFits)=2,msize(SappTEMrebuild_crop_WesslauFits)=2
	ModifyGraph msize(EET_ROIed_WesslauFits)=2,msize(thoSflow_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(METmetolu2_ROIed_WesslauFits)=2,msize(JHU102516_ROIed_WesslauFits)=9
	ModifyGraph msize(TL5b_ROIed_ROIed_ROI_WesslauFits)=2,msize(GRAmetolu_ROIed_WesslauFits)=2
	ModifyGraph msize(QUE_ROIed_WesslauFits)=2,msize(Orgueil2018_ROIed_Cr_WesslauFits)=2
	ModifyGraph msize(JHU122016_ROIed_WesslauFits)=2,msize(EETmetolu_ROIed_ROIe_WesslauFits)=2
	ModifyGraph msize(QUEmetolu2_ROIed_WesslauFits)=2,msize(QUEmetolu_ROIed_WesslauFits)=2
	ModifyGraph msize(Renazzo50_500_ROIed__WesslauFits)=2,msize(MaiaCONGR34_ROIed_WesslauFits)=2
	ModifyGraph msize(TL11i_ROIed_ROIed_cr_WesslauFits)=2,msize(JHU120616_ROIed_ROIe_WesslauFits)=2
	ModifyGraph msize(METmetolu_ROIed_WesslauFits)=2,msize(NMP20_ROIed_WesslauFits)=2
	ModifyGraph msize(Murch_Neg_UF_ROIed_N_WesslauFits)=2,msize(Fre211b_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(EETmetolu2_ROIed_WesslauFits)=2,msize(MaiaCO28_ROIed_WesslauFits)=2
	ModifyGraph msize(NMP5_ROIed_WesslauFits)=2,msize(thoDflow_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(GRAmetolu2_ROIed_WesslauFits)=2,msize(MACR5_ROIed_WesslauFits)=2
	ModifyGraph msize(JunkoIsa2_ROIed_ROIe_WesslauFits)=2,msize(NWApos_UF_ROIed_WesslauFits)=2
	ModifyGraph msize(MACR10_ROIed_WesslauFits)=2,msize(MaiaCOGR30_ROIed_WesslauFits)=2
	ModifyGraph msize(Murch_clean_ROIed_WesslauFits)=9,msize(MACR15_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(Fre311b_ROIed_WesslauFits)=2,msize(MACR20_ROIed_ROIed_WesslauFits)=2
	ModifyGraph msize(MVK5_ROIed_WesslauFits)=2,msize(Fre1011b_ROIed_WesslauFits)=2
	ModifyGraph msize(MVK0_ROIed_WesslauFits)=2,msize(Isono61WN_ROIed_WesslauFits)=2
	ModifyGraph msize(Nogoya50_500_ROIed_R_WesslauFits)=2,msize(fresneau_ROIed_WesslauFits)=2
	ModifyGraph msize(SucreHCL_ROIed_recal_WesslauFits)=6,msize(Fre315b_ROIed_WesslauFits)=2
	ModifyGraph msize(NMP0_ROIed_WesslauFits)=2,msize(Isono64WON_ROIed_WesslauFits)=2
	ModifyGraph msize(Fre211bHYDRO_ROIed_R_WesslauFits)=2,msize(MACR0_ROIed_WesslauFits)=2
	ModifyGraph msize(steack_ROIed_ROIed_WesslauFits)=2,msize(PSOCLDI_ROIed_WesslauFits)=2
	ModifyGraph msize(MZUVtemLEGIT_ROIed_R_WesslauFits)=2
	ModifyGraph opaque(JHU112216_ROIed_WesslauFits)=1,opaque(JHU112916_ROIed_WesslauFits)=1
	ModifyGraph opaque(JHU102516_ROIed_WesslauFits)=1,opaque(JHU122016_ROIed_WesslauFits)=1
	ModifyGraph opaque(JHU120616_ROIed_ROIe_WesslauFits)=1
	ModifyGraph hideTrace(Toxi2_crop_ROIed_ROI_WesslauFits)=1,hideTrace(MZUVirrLEGIT_ROIed_WesslauFits)=1
	ModifyGraph hideTrace(LS_ROIed_WesslauFits)=1,hideTrace(OrgueilLDI_WesslauFits)=1
	ModifyGraph hideTrace(GRA_ROIed_WesslauFits)=1,hideTrace(MET_ROIed_WesslauFits)=1
	ModifyGraph hideTrace(psocmeto_uf_cal_WesslauFits)=1,hideTrace(EET_ROIed_WesslauFits)=1
	ModifyGraph hideTrace(QUE_ROIed_WesslauFits)=1,hideTrace(TL11i_ROIed_ROIed_cr_WesslauFits)=1
	ModifyGraph hideTrace(NMP20_ROIed_WesslauFits)=1,hideTrace(Murch_Neg_UF_ROIed_N_WesslauFits)=1
	ModifyGraph hideTrace(NWApos_UF_ROIed_WesslauFits)=1,hideTrace(SucreHCL_ROIed_recal_WesslauFits)=1
	ModifyGraph hideTrace(steack_ROIed_ROIed_WesslauFits)=1,hideTrace(PSOCLDI_ROIed_WesslauFits)=1
	ModifyGraph hideTrace(MZUVtemLEGIT_ROIed_R_WesslauFits)=1
	ModifyGraph useMrkStrokeRGB=1
	ModifyGraph muloffset(Toxi2_ROIed_ROIed_WesslauFits)={1000000000,0},muloffset(Toxi3_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Toxi2_crop_ROIed_ROI_WesslauFits)={1000000000,0},muloffset(TOXI1_ROIed_ROIed_RO_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(MZUVirrLEGIT_ROIed_WesslauFits)={1000000000,0},muloffset(LS_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(OrgueilLDI_WesslauFits)={1000000000,0},muloffset(JHU112216_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Ivuna50_500_ROIed_RO_WesslauFits)={1000000000,0},muloffset(PicachuP072A_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(GRA_ROIed_WesslauFits)={1000000000,0},muloffset(PicachuP072B_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Sapp4_ROIed_ROIed_WesslauFits)={1000000000,0},muloffset(CBokk50_500_ROIed_RO_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(JHU112916_ROIed_WesslauFits)={1000000000,0},muloffset(hortholin_ROIed_ROIe_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Vass_ROIed_ROIed_ROI_WesslauFits)={1000000000,0},muloffset(MET_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(psocmeto_uf_cal_WesslauFits)={1000000000,0},muloffset(SappTEMrebuild_crop_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(EET_ROIed_WesslauFits)={1000000000,0},muloffset(thoSflow_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(METmetolu2_ROIed_WesslauFits)={1000000000,0},muloffset(JHU102516_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(TL5b_ROIed_ROIed_ROI_WesslauFits)={1000000000,0},muloffset(GRAmetolu_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(QUE_ROIed_WesslauFits)={1000000000,0},muloffset(Orgueil2018_ROIed_Cr_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(JHU122016_ROIed_WesslauFits)={1000000000,0},muloffset(EETmetolu_ROIed_ROIe_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(QUEmetolu2_ROIed_WesslauFits)={1000000000,0},muloffset(QUEmetolu_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Renazzo50_500_ROIed__WesslauFits)={1000000000,0},muloffset(MaiaCONGR34_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(TL11i_ROIed_ROIed_cr_WesslauFits)={1000000000,0},muloffset(JHU120616_ROIed_ROIe_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(METmetolu_ROIed_WesslauFits)={1000000000,0},muloffset(NMP20_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Murch_Neg_UF_ROIed_N_WesslauFits)={1000000000,0},muloffset(Fre211b_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(EETmetolu2_ROIed_WesslauFits)={1000000000,0},muloffset(MaiaCO28_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(NMP5_ROIed_WesslauFits)={1000000000,0},muloffset(thoDflow_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(GRAmetolu2_ROIed_WesslauFits)={1000000000,0},muloffset(MACR5_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(JunkoIsa2_ROIed_ROIe_WesslauFits)={1000000000,0},muloffset(NWApos_UF_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(MACR10_ROIed_WesslauFits)={1000000000,0},muloffset(MaiaCOGR30_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Murch_clean_ROIed_WesslauFits)={1000000000,0},muloffset(MACR15_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Fre311b_ROIed_WesslauFits)={1000000000,0},muloffset(MACR20_ROIed_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(MVK5_ROIed_WesslauFits)={1000000000,0},muloffset(Fre1011b_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(MVK0_ROIed_WesslauFits)={1000000000,0},muloffset(Isono61WN_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Nogoya50_500_ROIed_R_WesslauFits)={1000000000,0},muloffset(fresneau_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(SucreHCL_ROIed_recal_WesslauFits)={1000000000,0},muloffset(Fre315b_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(NMP0_ROIed_WesslauFits)={1000000000,0},muloffset(Isono64WON_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(Fre211bHYDRO_ROIed_R_WesslauFits)={1000000000,0},muloffset(MACR0_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(steack_ROIed_ROIed_WesslauFits)={1000000000,0},muloffset(PSOCLDI_ROIed_WesslauFits)={1000000000,0}
	ModifyGraph muloffset(MZUVtemLEGIT_ROIed_R_WesslauFits)={1000000000,0}
	ModifyGraph zmrkSize(Toxi2_ROIed_ROIed_WesslauFits)={Toxi2_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Toxi3_ROIed_ROIed_WesslauFits)={Toxi3_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Toxi2_crop_ROIed_ROI_WesslauFits)={Toxi2_crop_ROIed_ROI_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(TOXI1_ROIed_ROIed_RO_WesslauFits)={TOXI1_ROIed_ROIed_RO_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MZUVirrLEGIT_ROIed_WesslauFits)={MZUVirrLEGIT_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(LS_ROIed_WesslauFits)={LS_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(OrgueilLDI_WesslauFits)={OrgueilLDI_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JHU112216_ROIed_WesslauFits)={JHU112216_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Ivuna50_500_ROIed_RO_WesslauFits)={Ivuna50_500_ROIed_RO_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(PicachuP072A_ROIed_WesslauFits)={PicachuP072A_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(GRA_ROIed_WesslauFits)={GRA_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(PicachuP072B_ROIed_WesslauFits)={PicachuP072B_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Sapp4_ROIed_ROIed_WesslauFits)={Sapp4_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(CBokk50_500_ROIed_RO_WesslauFits)={CBokk50_500_ROIed_RO_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JHU112916_ROIed_WesslauFits)={JHU112916_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(hortholin_ROIed_ROIe_WesslauFits)={hortholin_ROIed_ROIe_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Vass_ROIed_ROIed_ROI_WesslauFits)={Vass_ROIed_ROIed_ROI_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MET_ROIed_WesslauFits)={MET_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(psocmeto_uf_cal_WesslauFits)={psocmeto_uf_cal_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(SappTEMrebuild_crop_WesslauFits)={SappTEMrebuild_crop_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(EET_ROIed_WesslauFits)={EET_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(thoSflow_ROIed_ROIed_WesslauFits)={thoSflow_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(METmetolu2_ROIed_WesslauFits)={METmetolu2_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JHU102516_ROIed_WesslauFits)={JHU102516_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(TL5b_ROIed_ROIed_ROI_WesslauFits)={TL5b_ROIed_ROIed_ROI_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(GRAmetolu_ROIed_WesslauFits)={GRAmetolu_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(QUE_ROIed_WesslauFits)={QUE_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Orgueil2018_ROIed_Cr_WesslauFits)={Orgueil2018_ROIed_Cr_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JHU122016_ROIed_WesslauFits)={JHU122016_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(EETmetolu_ROIed_ROIe_WesslauFits)={EETmetolu_ROIed_ROIe_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(QUEmetolu2_ROIed_WesslauFits)={QUEmetolu2_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(QUEmetolu_ROIed_WesslauFits)={QUEmetolu_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Renazzo50_500_ROIed__WesslauFits)={Renazzo50_500_ROIed__WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MaiaCONGR34_ROIed_WesslauFits)={MaiaCONGR34_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(TL11i_ROIed_ROIed_cr_WesslauFits)={TL11i_ROIed_ROIed_cr_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JHU120616_ROIed_ROIe_WesslauFits)={JHU120616_ROIed_ROIe_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(METmetolu_ROIed_WesslauFits)={METmetolu_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(NMP20_ROIed_WesslauFits)={NMP20_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Murch_Neg_UF_ROIed_N_WesslauFits)={Murch_Neg_UF_ROIed_N_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Fre211b_ROIed_ROIed_WesslauFits)={Fre211b_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(EETmetolu2_ROIed_WesslauFits)={EETmetolu2_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MaiaCO28_ROIed_WesslauFits)={MaiaCO28_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(NMP5_ROIed_WesslauFits)={NMP5_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(thoDflow_ROIed_ROIed_WesslauFits)={thoDflow_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(GRAmetolu2_ROIed_WesslauFits)={GRAmetolu2_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MACR5_ROIed_WesslauFits)={MACR5_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(JunkoIsa2_ROIed_ROIe_WesslauFits)={JunkoIsa2_ROIed_ROIe_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(NWApos_UF_ROIed_WesslauFits)={NWApos_UF_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MACR10_ROIed_WesslauFits)={MACR10_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MaiaCOGR30_ROIed_WesslauFits)={MaiaCOGR30_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Murch_clean_ROIed_WesslauFits)={Murch_clean_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MACR15_ROIed_ROIed_WesslauFits)={MACR15_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Fre311b_ROIed_WesslauFits)={Fre311b_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MACR20_ROIed_ROIed_WesslauFits)={MACR20_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MVK5_ROIed_WesslauFits)={MVK5_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Fre1011b_ROIed_WesslauFits)={Fre1011b_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MVK0_ROIed_WesslauFits)={MVK0_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Isono61WN_ROIed_WesslauFits)={Isono61WN_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Nogoya50_500_ROIed_R_WesslauFits)={Nogoya50_500_ROIed_R_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(fresneau_ROIed_WesslauFits)={fresneau_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(SucreHCL_ROIed_recal_WesslauFits)={SucreHCL_ROIed_recal_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Fre315b_ROIed_WesslauFits)={Fre315b_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(NMP0_ROIed_WesslauFits)={NMP0_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Isono64WON_ROIed_WesslauFits)={Isono64WON_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(Fre211bHYDRO_ROIed_R_WesslauFits)={Fre211bHYDRO_ROIed_R_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MACR0_ROIed_WesslauFits)={MACR0_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(steack_ROIed_ROIed_WesslauFits)={steack_ROIed_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(PSOCLDI_ROIed_WesslauFits)={PSOCLDI_ROIed_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zmrkSize(MZUVtemLEGIT_ROIed_R_WesslauFits)={MZUVtemLEGIT_ROIed_R_WesslauFits[13][*],-4,-1,6,1}
	ModifyGraph zColor(Toxi2_ROIed_ROIed_WesslauFits)={Toxi2_ROIed_ROIed_WesslauFits[19][*],-1.45,0.5,dBZ21}
	ModifyGraph zColor(Toxi3_ROIed_ROIed_WesslauFits)={Toxi3_ROIed_ROIed_WesslauFits[19][*],-1.45,0.5,dBZ21}
	ModifyGraph zColor(Toxi2_crop_ROIed_ROI_WesslauFits)={Toxi2_crop_ROIed_ROI_WesslauFits[19][*],-1.45,0.5,dBZ21}
	ModifyGraph zColor(TOXI1_ROIed_ROIed_RO_WesslauFits)={TOXI1_ROIed_ROIed_RO_WesslauFits[19][*],-1.45,0.5,dBZ21}
	ModifyGraph zColor(MZUVirrLEGIT_ROIed_WesslauFits)={MZUVirrLEGIT_ROIed_WesslauFits[19][*],-1.5,0.5,Rainbow16}
	ModifyGraph zColor(LS_ROIed_WesslauFits)={LS_ROIed_WesslauFits[19][*],0,0.5,Spectrum}
	ModifyGraph zColor(OrgueilLDI_WesslauFits)={OrgueilLDI_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(JHU112216_ROIed_WesslauFits)={JHU112216_ROIed_WesslauFits[19][*],-0.1,0.5,Grays16}
	ModifyGraph zColor(Ivuna50_500_ROIed_RO_WesslauFits)={Ivuna50_500_ROIed_RO_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(PicachuP072A_ROIed_WesslauFits)={PicachuP072A_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(GRA_ROIed_WesslauFits)={GRA_ROIed_WesslauFits[19][*],-0.25,0.5,Rainbow}
	ModifyGraph zColor(PicachuP072B_ROIed_WesslauFits)={PicachuP072B_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(Sapp4_ROIed_ROIed_WesslauFits)={Sapp4_ROIed_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(CBokk50_500_ROIed_RO_WesslauFits)={CBokk50_500_ROIed_RO_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(JHU112916_ROIed_WesslauFits)={JHU112916_ROIed_WesslauFits[19][*],-0.1,0.5,Grays16}
	ModifyGraph zColor(hortholin_ROIed_ROIe_WesslauFits)={hortholin_ROIed_ROIe_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(Vass_ROIed_ROIed_ROI_WesslauFits)={Vass_ROIed_ROIed_ROI_WesslauFits[19][*],0,0.5,RedWhiteGreen}
	ModifyGraph zColor(MET_ROIed_WesslauFits)={MET_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(psocmeto_uf_cal_WesslauFits)={psocmeto_uf_cal_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(SappTEMrebuild_crop_WesslauFits)={SappTEMrebuild_crop_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(EET_ROIed_WesslauFits)={EET_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(thoSflow_ROIed_ROIed_WesslauFits)={thoSflow_ROIed_ROIed_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(METmetolu2_ROIed_WesslauFits)={METmetolu2_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(JHU102516_ROIed_WesslauFits)={JHU102516_ROIed_WesslauFits[19][*],-0.1,0.5,Grays16}
	ModifyGraph zColor(TL5b_ROIed_ROIed_ROI_WesslauFits)={TL5b_ROIed_ROIed_ROI_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(GRAmetolu_ROIed_WesslauFits)={GRAmetolu_ROIed_WesslauFits[19][*],-0.25,0.5,Grays}
	ModifyGraph zColor(QUE_ROIed_WesslauFits)={QUE_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(Orgueil2018_ROIed_Cr_WesslauFits)={Orgueil2018_ROIed_Cr_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(JHU122016_ROIed_WesslauFits)={JHU122016_ROIed_WesslauFits[19][*],-0.1,0.5,Grays16}
	ModifyGraph zColor(EETmetolu_ROIed_ROIe_WesslauFits)={EETmetolu_ROIed_ROIe_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(QUEmetolu2_ROIed_WesslauFits)={QUEmetolu2_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(QUEmetolu_ROIed_WesslauFits)={QUEmetolu_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(Renazzo50_500_ROIed__WesslauFits)={Renazzo50_500_ROIed__WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(MaiaCONGR34_ROIed_WesslauFits)={MaiaCONGR34_ROIed_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(TL11i_ROIed_ROIed_cr_WesslauFits)={TL11i_ROIed_ROIed_cr_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(JHU120616_ROIed_ROIe_WesslauFits)={JHU120616_ROIed_ROIe_WesslauFits[19][*],-0.1,0.5,Grays16}
	ModifyGraph zColor(METmetolu_ROIed_WesslauFits)={METmetolu_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(NMP20_ROIed_WesslauFits)={NMP20_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(Murch_Neg_UF_ROIed_N_WesslauFits)={Murch_Neg_UF_ROIed_N_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(Fre211b_ROIed_ROIed_WesslauFits)={Fre211b_ROIed_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(EETmetolu2_ROIed_WesslauFits)={EETmetolu2_ROIed_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(MaiaCO28_ROIed_WesslauFits)={MaiaCO28_ROIed_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(NMP5_ROIed_WesslauFits)={NMP5_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(thoDflow_ROIed_ROIed_WesslauFits)={thoDflow_ROIed_ROIed_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(GRAmetolu2_ROIed_WesslauFits)={GRAmetolu2_ROIed_WesslauFits[19][*],-0.25,0.5,Grays}
	ModifyGraph zColor(MACR5_ROIed_WesslauFits)={MACR5_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(JunkoIsa2_ROIed_ROIe_WesslauFits)={JunkoIsa2_ROIed_ROIe_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(NWApos_UF_ROIed_WesslauFits)={NWApos_UF_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(MACR10_ROIed_WesslauFits)={MACR10_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(MaiaCOGR30_ROIed_WesslauFits)={MaiaCOGR30_ROIed_WesslauFits[19][*],-0.5,0.5,Rainbow}
	ModifyGraph zColor(Murch_clean_ROIed_WesslauFits)={Murch_clean_ROIed_WesslauFits[19][*],-1.5,0.5,Rainbow16}
	ModifyGraph zColor(MACR15_ROIed_ROIed_WesslauFits)={MACR15_ROIed_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(Fre311b_ROIed_WesslauFits)={Fre311b_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(MACR20_ROIed_ROIed_WesslauFits)={MACR20_ROIed_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(MVK5_ROIed_WesslauFits)={MVK5_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(Fre1011b_ROIed_WesslauFits)={Fre1011b_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(MVK0_ROIed_WesslauFits)={MVK0_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(Isono61WN_ROIed_WesslauFits)={Isono61WN_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(Nogoya50_500_ROIed_R_WesslauFits)={Nogoya50_500_ROIed_R_WesslauFits[19][*],0,0.5,Grays}
	ModifyGraph zColor(fresneau_ROIed_WesslauFits)={fresneau_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(SucreHCL_ROIed_recal_WesslauFits)={SucreHCL_ROIed_recal_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(Fre315b_ROIed_WesslauFits)={Fre315b_ROIed_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(NMP0_ROIed_WesslauFits)={NMP0_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(Isono64WON_ROIed_WesslauFits)={Isono64WON_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(Fre211bHYDRO_ROIed_R_WesslauFits)={Fre211bHYDRO_ROIed_R_WesslauFits[19][*],0,0.5,CyanMagenta}
	ModifyGraph zColor(MACR0_ROIed_WesslauFits)={MACR0_ROIed_WesslauFits[19][*],0,0.5,RedWhiteBlue}
	ModifyGraph zColor(steack_ROIed_ROIed_WesslauFits)={steack_ROIed_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(PSOCLDI_ROIed_WesslauFits)={PSOCLDI_ROIed_WesslauFits[19][*],0,0.5,BlueRedGreen}
	ModifyGraph zColor(MZUVtemLEGIT_ROIed_R_WesslauFits)={MZUVtemLEGIT_ROIed_R_WesslauFits[19][*],-1.5,0.5,Rainbow16}
	ModifyGraph zColorMax(Toxi2_ROIed_ROIed_WesslauFits)=NaN,zColorMax(Toxi3_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Toxi2_crop_ROIed_ROI_WesslauFits)=NaN,zColorMax(TOXI1_ROIed_ROIed_RO_WesslauFits)=NaN
	ModifyGraph zColorMax(MZUVirrLEGIT_ROIed_WesslauFits)=NaN,zColorMax(LS_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(OrgueilLDI_WesslauFits)=NaN,zColorMax(JHU112216_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Ivuna50_500_ROIed_RO_WesslauFits)=NaN,zColorMax(PicachuP072A_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(GRA_ROIed_WesslauFits)=NaN,zColorMax(PicachuP072B_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Sapp4_ROIed_ROIed_WesslauFits)=NaN,zColorMax(CBokk50_500_ROIed_RO_WesslauFits)=NaN
	ModifyGraph zColorMax(JHU112916_ROIed_WesslauFits)=NaN,zColorMax(hortholin_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMax(Vass_ROIed_ROIed_ROI_WesslauFits)=NaN,zColorMax(MET_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(psocmeto_uf_cal_WesslauFits)=NaN,zColorMax(SappTEMrebuild_crop_WesslauFits)=NaN
	ModifyGraph zColorMax(EET_ROIed_WesslauFits)=NaN,zColorMax(thoSflow_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(METmetolu2_ROIed_WesslauFits)=NaN,zColorMax(JHU102516_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(TL5b_ROIed_ROIed_ROI_WesslauFits)=NaN,zColorMax(GRAmetolu_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(QUE_ROIed_WesslauFits)=NaN,zColorMax(Orgueil2018_ROIed_Cr_WesslauFits)=NaN
	ModifyGraph zColorMax(JHU122016_ROIed_WesslauFits)=NaN,zColorMax(EETmetolu_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMax(QUEmetolu2_ROIed_WesslauFits)=NaN,zColorMax(QUEmetolu_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Renazzo50_500_ROIed__WesslauFits)=NaN,zColorMax(MaiaCONGR34_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(TL11i_ROIed_ROIed_cr_WesslauFits)=NaN,zColorMax(JHU120616_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMax(METmetolu_ROIed_WesslauFits)=NaN,zColorMax(NMP20_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Murch_Neg_UF_ROIed_N_WesslauFits)=NaN,zColorMax(Fre211b_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(EETmetolu2_ROIed_WesslauFits)=NaN,zColorMax(MaiaCO28_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(NMP5_ROIed_WesslauFits)=NaN,zColorMax(thoDflow_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(GRAmetolu2_ROIed_WesslauFits)=NaN,zColorMax(MACR5_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(JunkoIsa2_ROIed_ROIe_WesslauFits)=NaN,zColorMax(NWApos_UF_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(MACR10_ROIed_WesslauFits)=NaN,zColorMax(MaiaCOGR30_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Murch_clean_ROIed_WesslauFits)=NaN,zColorMax(MACR15_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Fre311b_ROIed_WesslauFits)=NaN,zColorMax(MACR20_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(MVK5_ROIed_WesslauFits)=NaN,zColorMax(Fre1011b_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(MVK0_ROIed_WesslauFits)=NaN,zColorMax(Isono61WN_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Nogoya50_500_ROIed_R_WesslauFits)=NaN,zColorMax(fresneau_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(SucreHCL_ROIed_recal_WesslauFits)=NaN,zColorMax(Fre315b_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(NMP0_ROIed_WesslauFits)=NaN,zColorMax(Isono64WON_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(Fre211bHYDRO_ROIed_R_WesslauFits)=NaN,zColorMax(MACR0_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(steack_ROIed_ROIed_WesslauFits)=NaN,zColorMax(PSOCLDI_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMax(MZUVtemLEGIT_ROIed_R_WesslauFits)=NaN
	ModifyGraph zColorMin(Toxi2_ROIed_ROIed_WesslauFits)=NaN,zColorMin(Toxi3_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Toxi2_crop_ROIed_ROI_WesslauFits)=NaN,zColorMin(TOXI1_ROIed_ROIed_RO_WesslauFits)=NaN
	ModifyGraph zColorMin(MZUVirrLEGIT_ROIed_WesslauFits)=NaN,zColorMin(LS_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(OrgueilLDI_WesslauFits)=NaN,zColorMin(JHU112216_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Ivuna50_500_ROIed_RO_WesslauFits)=NaN,zColorMin(PicachuP072A_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(GRA_ROIed_WesslauFits)=NaN,zColorMin(PicachuP072B_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Sapp4_ROIed_ROIed_WesslauFits)=NaN,zColorMin(CBokk50_500_ROIed_RO_WesslauFits)=NaN
	ModifyGraph zColorMin(JHU112916_ROIed_WesslauFits)=NaN,zColorMin(hortholin_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMin(Vass_ROIed_ROIed_ROI_WesslauFits)=NaN,zColorMin(MET_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(psocmeto_uf_cal_WesslauFits)=NaN,zColorMin(SappTEMrebuild_crop_WesslauFits)=NaN
	ModifyGraph zColorMin(EET_ROIed_WesslauFits)=NaN,zColorMin(thoSflow_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(METmetolu2_ROIed_WesslauFits)=NaN,zColorMin(JHU102516_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(TL5b_ROIed_ROIed_ROI_WesslauFits)=NaN,zColorMin(GRAmetolu_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(QUE_ROIed_WesslauFits)=NaN,zColorMin(Orgueil2018_ROIed_Cr_WesslauFits)=NaN
	ModifyGraph zColorMin(JHU122016_ROIed_WesslauFits)=NaN,zColorMin(EETmetolu_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMin(QUEmetolu2_ROIed_WesslauFits)=NaN,zColorMin(QUEmetolu_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Renazzo50_500_ROIed__WesslauFits)=NaN,zColorMin(MaiaCONGR34_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(TL11i_ROIed_ROIed_cr_WesslauFits)=NaN,zColorMin(JHU120616_ROIed_ROIe_WesslauFits)=NaN
	ModifyGraph zColorMin(METmetolu_ROIed_WesslauFits)=NaN,zColorMin(NMP20_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Murch_Neg_UF_ROIed_N_WesslauFits)=NaN,zColorMin(Fre211b_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(EETmetolu2_ROIed_WesslauFits)=NaN,zColorMin(MaiaCO28_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(NMP5_ROIed_WesslauFits)=NaN,zColorMin(thoDflow_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(GRAmetolu2_ROIed_WesslauFits)=NaN,zColorMin(MACR5_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(JunkoIsa2_ROIed_ROIe_WesslauFits)=NaN,zColorMin(NWApos_UF_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(MACR10_ROIed_WesslauFits)=NaN,zColorMin(MaiaCOGR30_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Murch_clean_ROIed_WesslauFits)=NaN,zColorMin(MACR15_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Fre311b_ROIed_WesslauFits)=NaN,zColorMin(MACR20_ROIed_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(MVK5_ROIed_WesslauFits)=NaN,zColorMin(Fre1011b_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(MVK0_ROIed_WesslauFits)=NaN,zColorMin(Isono61WN_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Nogoya50_500_ROIed_R_WesslauFits)=NaN,zColorMin(fresneau_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(SucreHCL_ROIed_recal_WesslauFits)=NaN,zColorMin(Fre315b_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(NMP0_ROIed_WesslauFits)=NaN,zColorMin(Isono64WON_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(Fre211bHYDRO_ROIed_R_WesslauFits)=NaN,zColorMin(MACR0_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(steack_ROIed_ROIed_WesslauFits)=NaN,zColorMin(PSOCLDI_ROIed_WesslauFits)=NaN
	ModifyGraph zColorMin(MZUVtemLEGIT_ROIed_R_WesslauFits)=NaN
	ModifyGraph log=1
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=4
	ModifyGraph axThick=2
	Label left "Growth factor (arbitrary unit)"
	Label bottom "Average number of reactions per molecule"
	SetAxis left 1.00070091801331,1.11398740789127
	SetAxis bottom 3,407.97593
	Legend/C/N=text0_1/J/A=MC/X=11.64/Y=44.18 "\\Z20\\s(Toxi2_ROIed_ROIed_WesslauFits) Water-insoluble Murchison"
	AppendText "\\s(Toxi2_ROIed_ROIed_WesslauFits) Carbonaceous chondrites extracts\r\t\tMethanol and Methanol/Toluene extracts"
	Legend/C/N=text0_2/J/A=MC/X=24.19/Y=27.86 "\\Z20\\s(Toxi2_ROIed_ROIed_WesslauFits) Gas phase nebular simulants"
	AppendText "\t\tHigh hydrogen content (50% < f(H\\B2\\M\\Z20) < 62%) ~1800K\r\t\t\\f02Nebulotron\\f00 experiments c.f. Bekaert et alii, 2018"
	AppendText "\\s(Toxi2_ROIed_ROIed_WesslauFits) Gas phase exoplanet simulants\r\t\tHigh hydrogen content (x(H) 100 to 1000 times solar) <600K"
	AppendText "\t\t\\f02PHAZER\\f00 experiments c.f. Horst et alii, 2018"
	Legend/C/N=text0_3/J/A=MC/X=26.61/Y=10.51 "\\Z20\\s(Toxi2_ROIed_ROIed_WesslauFits) Titan atmosphere and dry Solar Nebula simulants"
	AppendText "\t\tLow hydrogen content (f(H\\B2\\M\\Z20)=0, x(H)<<50%)\r\t\tc.f. Kuga et alii 2014, 2015, 2017 and Horst et alii 2013"
	Legend/C/N=text0_4/J/A=MC/X=32.58/Y=-11.82 "\\Z20\\s(Toxi2_ROIed_ROIed_WesslauFits) Aqueous chemistry\r\t\tc.f. Isono et alii, 2019, Renard et alii, 2013"
	AppendText "\t\t\tVinogradoff et alii, 2018"
	Legend/C/N=text0_5/J/A=MC/X=30.01/Y=-1.03 "\\Z20\\s(Toxi2_ROIed_ROIed_WesslauFits) Irradiated ices, cometary simulants"
	AppendText "\t\tc.f. Fresneau et alii, 2017 and Piani et alii, 2017"
EndMacro

function/WAVE SetUnionProx(w)
wave w
variable n=dimsize(w,0),k
Duplicate/O/FREE w,parent
parent=x
parent[0]=1
parent[n-1]=n-2
for(k=1;k<(n-1);k+=1)
	if(w[k]>(w[k+1]+w[k-1])/2)
		parent[k]=find(k+1,parent)
	else
		parent[k]=find(k-1,parent)
	endif
endfor
return parent
end

function/WAVE differn(w,n)
wave w
variable n
variable m=dimsize(w,0)
make/O/D/N=(m,n) diffsucc=1
diffsucc[1,m-2][0]=(w[x+1]-w[x-1])/2
diffsucc[0][0]=w[1]-w[0]
diffsucc[m-1][0]=w[m-1]-w[m-2]
variable k,j
for(k=1;k<n;k+=1)
	for(j=k;j<(m-k);j+=1)
		diffsucc[j][k]=diffsucc[j+1][k-1]-diffsucc[j-1][k-1]
	endfor
endfor
return diffsucc
end

function sixChroma2roi(name)
string name
name=name+"_6chroma"
wave sixChroma=$name
variable n= dimsize(sixChroma,0)
make/O/D/N=(n) roi0,roi1, roi2, roi3
roi0=sixChroma[x][6]
roi1=sixChroma[x][1]
roi2=defmass(roi0)
roi3=sixChroma[x][7]
sort roi0,roi0,roi1,roi2,roi3
end

function roi2active()
wave roi0,roi1,roi2,roi3
Duplicate/O roi0 wave0stripped
Duplicate/O roi1 wave1stripped
Duplicate/O roi2 wave2stripped
Duplicate/O roi3 wave3stripped
end

function speClu4chromato(nb)
variable nb
wave roi3,roi0
variable n=numpnts(roi0)
make/O/D/N=(n,2)/FREE obs
obs[][0]=0//roi0[x]
obs[][1]=roi3[x]
Make/O/D/N=1/FREE grPar=pi
spectralClu(obs,"RNG",grPar,2,max(nb,2))
wave W_KMMembers,M_eigenVectors
if(nb==0)
	wave Laplace
	variable ncla=optiHyperParamSpeclu(obs,Laplace)
	print ncla
	classify(M_eigenVectors,2,ncla)
	obs[][0]=roi0[x]
	silhouette(obs,W_KMMembers)
	wave classesChiefs,classesEffectifs
	Duplicate/O W_KMMembers Clu_Cla3
	Duplicate/O classesChiefs Clu_Cen3
	Duplicate/O classesEffectifs Clu_Car3
	Killwaves classesChiefs,classesEffectifs,Laplace
else
	obs[][0]=roi0[x]
	silhouette(obs,W_KMMembers)
endif
end

function kmean4chromato(nb)
variable nb
wave roi3,roi0
variable n=numpnts(roi0)
make/O/D/N=(n,2)/FREE obs
obs[][0]=0//roi0[x]
obs[][1]=roi3[x]
if(nb)
else
	nb=optiHyperParamKmean(obs)
	print nb
endif
classify(obs,2,nb)
wave W_KMMembers
obs[][0]=roi0[x]
print silhouette(obs,W_KMMembers)
wave classesChiefs,classesEffectifs
Duplicate/O W_KMMembers Clu_Cla3
Duplicate/O classesChiefs Clu_Cen3
Duplicate/O classesEffectifs Clu_Car3
Killwaves classesChiefs,classesEffectifs
end

function optiHyperParamKmean(obs)
wave obs
variable i,n,dim
n=dimsize(obs,0)
dim=dimsize(obs,1)
Make/O/D/N=(n) silres=-1
for(i=2;i<n;i+=1)
	classify(obs,dim,i)
	wave W_KMMembers
	silres[i]=silhouette(obs,W_KMMembers)
endfor
variable xeu=wavemax(silres)
extract/O/FREE/INDX silres,xeuz,silres==xeu
return wavemin(xeuz)
end

function spectralClu(obs,graphType,graphParameters,kvec,ncla)
wave obs,graphParameters
variable kvec,ncla
string graphtype
Duplicate/O/FREE initiateAdja(Obs) Adja
Adja=1
Duplicate/O/FREE calcDistances(obs,Adja) distances
strswitch(graphType)	// string switch
	case "RNG":
		Duplicate/O/FREE Relativeneighborhood(obs,distances) Adja
		break
	case "betaSkeleton":
		Duplicate/O/FREE BetaSkeleton(Obs,graphParameters[0]) Adja
		break
	default:
		Adja=1
		print "default case"
endswitch
Adja*=exp(-distances)
Duplicate/O calcLaplace(Adja) Laplace
MatrixEigenV/SYM/EVEC Laplace
wave M_eigenVectors
classify(M_eigenVectors,kvec,ncla)
wave W_KMMembers
Duplicate/O W_KMMembers Clu_Cla3
silhouette(obs,W_KMMembers)
wave classesChiefs,classesEffectifs
Duplicate/O classesChiefs Clu_Cen3
Duplicate/O classesEffectifs Clu_Car3
Killwaves classesChiefs,classesEffectifs
end

function silhouette(obs,classLab)
wave obs,classLab
variable n=dimsize(obs,0)
variable dim=dimsize(obs,1)
variable nb=wavemax(classLab)+1
variable k,i,j,l,nk,ni
Make/O/D/N=(nb,dim) classesChiefs=0
Make/O/D/N=(nb) classesEffectifs=0
Make/O/D/N=(nb)/FREE closestChief
Make/O/D/N=(n)/FREE meanIntraDis=0, meanSecDis=0,silhou
variable dis=0, ref=inf
for(k=0;k<nb;k+=1)
	extract/O/FREE/INDX classLab, kindex, classLab==k
	nk=numpnts(kindex)
	classesEffectifs[k]=nk
	for(i=0;i<nk;i+=1)
		classesChiefs[k][]+=obs[kindex[i]][y]/nk
		for(j=0;j<nk;j+=1)
			dis=0
			for(l=0;l<dim;l+=1)
				dis+=(obs[kindex[i]][l]-obs[kindex[j]][l])^2
			endfor
			dis=sqrt(dis)
			meanIntraDis[kindex[i]]+=dis/(nk)
		endfor
	endfor
endfor
for(k=0;k<nb;k+=1)
	ref=inf
	for(i=0;i<nb;i+=1)
		if(k!=i)
			dis=0
			for(l=0;l<dim;l+=1)
				dis+=(classesChiefs[k][l]-classesChiefs[i][l])^2
			endfor
			dis=sqrt(dis)
			if(dis<ref)
				closestChief[k]=i
				ref=dis
			endif
		endif
	endfor
endfor
for(k=0;k<nb;k+=1)
	extract/O/FREE/INDX classLab, kindex, classLab==k
	extract/O/FREE/INDX classLab, iindex, classLab==closestChief[k]
	nk=classesEffectifs[k]
	ni=classesEffectifs[closestChief[k]]
	for(i=0;i<nk;i+=1)
		for(j=0;j<ni;j+=1)
			dis=0
			for(l=0;l<dim;l+=1)
				dis+=(obs[kindex[i]][l]-obs[iindex[j]][l])^2
			endfor
			dis=sqrt(dis)
			meanSecDis[kindex[i]]+=dis/(ni+1)
		endfor
	endfor
endfor
silhou=(meanSecDis[x]-meanIntraDis[x])/max(meanSecDis[x],meanIntraDis[x])
return mean(silhou)
end

Function classify(eigen,kvec,ncla)
wave eigen
variable kvec,ncla
variable n=dimsize(eigen,0)
make/O/D/N=(kvec,n)/FREE pop
pop=eigen[y][x][0]
make/O/D/N=(kvec,ncla)/FREE iniclas
KMeans/INW=iniclas/OUT=2 pop
end

function/WAVE calcDistances(Obs,Adja)
wave obs,Adja
variable n=dimsize(obs,0)
variable dim=dimsize(obs,1)
make/O/D/N=(n,n)/FREE distances=Inf
variable k,j,contdim
for(k=0;k<n;k+=1)
	for(j=0;j<(k+1);j+=1)
		if(Adja[k][j])
			distances[k][j]=0
			for(contdim=0;contdim<dim;contdim+=1)
				distances[k][j]+=(obs[k][contdim]-obs[j][contdim])^2
			endfor
			distances[k][j]=sqrt(distances[k][j])
			distances[j][k]=distances[k][j]
		endif
	endfor
endfor
return Distances
end

function/wave Relativeneighborhood(obs,poids)//returns the adjacence matrix of the RNGraph for the poids metric in Obs set
wave obs,poids
variable n=dimsize(obs,0)
variable dim=dimsize(obs,1)
Make/O/D/N=(n,n)/FREE Adja=1
variable k,j,m,i,ref, dis1, dis2
//pour tous points
for(j=1;j<n;j+=1)//pour toute ligne de poids
	for(k=0;k<j;k+=1)//et toute colonne de poids sans doublons
		ref=poids[j][k]
		for(m=0;m<n;m+=1)
			dis1=poids[m][k]
			dis2=poids[j][m]
			if(dis1<ref && dis2<ref)
				if(dis1>0 && dis2>0)
					m=n
					Adja[j][k]=0
					Adja[k][j]=0
				endif
			endif
		endfor
	endfor
endfor
return Adja
end

function/WAVE BetaSkeleton(Obs,bet)
wave obs
variable bet
variable n=dimsize(obs,0)
variable dim=dimsize(obs,1)
Make/O/D/N=(n,n)/FREE Adja=1
Make/O/D/N=(dim)/FREE vec1,vec2,prod
variable k,j,m
variable angle
for(j=1;j<n;j+=1)//pour toute ligne de poids
	for(k=0;k<j;k+=1)//et toute colonne de poids sans doublons
		for(m=0;m<n;m+=1)
			vec1=obs[j][x]-obs[m][x]
			vec2=obs[k][x]-obs[m][x]
			prod=vec1*vec2
			angle=acos(sum(prod)/(norm(vec1)*norm(vec2)))
			if(angle>bet && (m!=k && m!=j))
				m=n
				Adja[j][k]=0
				Adja[k][j]=0
			endif
		endfor
	endfor
endfor
return Adja
end

function/wave calcLaplace(Adja)
wave adja
matrixop/O/FREE Degree=sumrows(Adja)
matrixOP/O/FREE laplace=diagonal(degree)
laplace-=adja
Redimension/D laplace
variable dd=norm(laplace)
laplace/=dd
return laplace
end

function optiHyperParamSpeclu(obs,laplace)
wave laplace,obs
variable i,n,dim
n=dimsize(obs,0)
dim=dimsize(obs,1)
Make/O/D/N=(n) silres=-1
MatrixEigenV/SYM/EVEC laplace
for(i=2;i<n;i+=1)
	wave M_eigenVectors
	classify(M_eigenVectors,dim,i)
	wave W_KMMembers
	silres[i]=silhouette(obs,W_KMMembers)
endfor
variable xeu=wavemax(silres)
extract/O/FREE/INDX silres,xeuz,silres==xeu
return wavemin(xeuz)
end

function fast3Dplot(name)
string name
name=name+"_6chroma"
wave sixChroma=$name
variable n=dimsize(sixChroma,0)
MAKE/O/D/N=(n,3) Scat3D,Int3D
Int3D[][0]=sixChroma[x][6]
Int3D[][1]=sixChroma[x][7]
Int3D[][2]=log(sixChroma[x][1])
Scat3D[][0]=sixChroma[x][6]
Scat3D[][1]=sixChroma[x][7]
Scat3D[][2]=defmass(sixChroma[x][6])//log(sixChroma[x][1])
end

function fast3Droi()
wave roi0,roi1,roi2,roi3
variable n=dimsize(roi0,0)
MAKE/O/D/N=(n,3) RoiInt3D, RoiDefmass3D
RoiInt3D[][0]=roi0[x]
RoiInt3D[][1]=roi3[x]
RoiInt3D[][2]=roi1[x]
RoiDefmass3D[][0]=roi0[x]
RoiDefmass3D[][1]=roi3[x]
RoiDefmass3D[][2]=roi2[x]
end

function splitCalibration(bagName,classLab,polyN)
string bagname
wave classLab
variable polyN
bagname="Ciblebag_"+bagname
wave Ciblebag=$bagname
//produire les deux fonctions de calibration
extract/O/FREE/INDX classLab,les0,classLab==0
extract/O/FREE/INDX classLab,les1,classLab==1
Make/O/D/N=(numpnts(les0)) lesy0,lesx0
Make/O/D/N=(numpnts(les1)) lesy1,lesx1
lesy0=ciblebag[les0[x]][1]-target2theo(ciblebag[les0[x]][1],ciblebag[les0[x]][3])
lesx0=ciblebag[les0[x]][1]
lesy1=ciblebag[les1[x]][1]-target2theo(ciblebag[les1[x]][1],ciblebag[les1[x]][3])
lesx1=ciblebag[les1[x]][1]
CurveFit/Q/M=2/W=0 poly (polyN+1), lesy0/X=lesx0/D
wave W_coef
Duplicate/O/FREE W_Coef FitCoef0
CurveFit/Q/M=2/W=0 poly (polyN+1), lesy1/X=lesx1/D
Duplicate/O/FREE W_Coef FitCoef1
Duplicate/O lesy0 fit0
Duplicate/O lesy1 fit1
fit0=poly(FitCoef0,lesx0)
fit1=poly(FitCoef1,lesx1)
display/K=1
appendtograph lesy0 vs lesx0
appendtograph lesy1 vs lesx1
appendtograph fit0 vs lesx0
appendtograph fit1 vs lesx1
modifygraph rgb(lesy1)=(0,0,65535)
modifygraph rgb(fit0)=(65535,0,52428)
modifygraph rgb(fit1)=(65535,0,52428)
//produire la classification par regression logistique
extract/O/FREE Ciblebag,lesmass,y==1
CurveFit/M=2/W=0 Sigmoid, classLab/X=lesmass/D
Duplicate/O W_Coef SigCoef
//corriger la masse de la donnée active
wave wave0stripped
Duplicate/O wave0stripped error
error=(1-(SigCoef[0] + SigCoef[1]/(1+exp(-(wave0stripped[x]-SigCoef[2])/SigCoef[3]))))*poly(FitCoef0,wave0stripped[x])+(SigCoef[0] + SigCoef[1]/(1+exp(-(wave0stripped[x]-SigCoef[2])/SigCoef[3])))*poly(FitCoef1,wave0stripped[x])
appendtograph error vs wave0stripped
modifygraph rgb(error)=(2,39321,1)
end


Window Plot3D_int() : GizmoPlot
	PauseUpdate; Silent 1		// building window...
	// Building Gizmo 8 window...
	NewGizmo/W=(30.75,56.75,709.5,689)
	ModifyGizmo startRecMacro=700
	ModifyGizmo scalingOption=63
	AppendToGizmo Scatter=root:Int3D,name=scatter0
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,3}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,0.1}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ markerCTab,dBZ14}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ CTABScaling,512}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ DropLines,16}
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisColor,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,ticks,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,ticks,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,ticks,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelText,"Mass"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelText,"Intensity"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelText,"Time"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,labelBillboarding,1}
	ModifyGizmo modifyObject=axes0,objectType=Axes,property={-1,Clipped,0}
	AppendToGizmo Scatter=root:RoiInt3D,name=scatter1
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ scatterColorType,0}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ size,0.3}
	ModifyGizmo ModifyObject=scatter1,objectType=scatter,property={ color,1,0,0,1}
	ModifyGizmo setDisplayList=0, object=scatter0
	ModifyGizmo setDisplayList=1, object=axes0
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ShowTools
	ModifyGizmo showInfo
	ModifyGizmo infoWindow={341,0,1166,284}
	ModifyGizmo endRecMacro
	ModifyGizmo SETQUATERNION={0.707107,-0.000000,-0.000000,0.707107}
EndMacro

Window Plot3D_DM() : GizmoPlot
	PauseUpdate; Silent 1		// building window...
	// Building Gizmo 8 window...
	NewGizmo/W=(914.25,530.75,1248.75,858.5)
	ModifyGizmo startRecMacro=700
	ModifyGizmo scalingOption=63
	AppendToGizmo Scatter=root:Scat3D,name=scatter0
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,3}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,0.15}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ markerCTab,ColdWarm}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ CTABScaling,66}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ MaxRGBA,1534,2.36943e-38,2.36943e-38,2.36943e-38,1}
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,gridType,41}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,gridPlaneColor,0.933333,0.933333,0.933333,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,gridLinesColor,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,gridPrimaryCount,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,gridSecondaryCount,5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisColor,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,ticks,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelText,"Mass"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelText,"Mass defect"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelText,"Time"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={5,visible,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={7,visible,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={11,visible,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,labelBillboarding,1}
	ModifyGizmo modifyObject=axes0,objectType=Axes,property={-1,Clipped,0}
	AppendToGizmo Path=root:Traceur_Chroma,name=path0
	ModifyGizmo ModifyObject=path0,objectType=path,property={ pathColorType,1}
	ModifyGizmo ModifyObject=path0,objectType=path,property={ lineWidthType,0}
	ModifyGizmo ModifyObject=path0,objectType=path,property={ pathColor,0,0,0,1}
	ModifyGizmo setDisplayList=0, object=scatter0
	ModifyGizmo setDisplayList=1, object=axes0
	ModifyGizmo setDisplayList=2, object=path0
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ShowTools
	ModifyGizmo showInfo
	ModifyGizmo infoWindow={929,1152,1530,1426}
	ModifyGizmo endRecMacro
	ModifyGizmo SETQUATERNION={0.369109,0.120770,0.294584,0.873153}
EndMacro


function ReformateAdja4chromato(adja,parent,roipnts,classif,cle)
wave adja,parent,roipnts,classif
variable cle
wave wave1stripped, wave0stripped,wave2stripped
//couper les ponts
variable n=dimsize(adja,0),k
variable m=numpnts(roipnts),i
for(i=0;i<m;i+=1)
	parent[roipnts[i]]=roipnts[i]
	for(k=0;k<n;k+=1)
		adja[roipnts[i]][k]=0
		adja[k][roipnts[i]]=0
	endfor
endfor
//reformer les chaines suivant la classif
variable nbcla=wavemax(classif)+1,j
for(j=0;j<nbcla;j+=1)
	extract/O/FREE roipnts,croproi,classif==j
	m=numpnts(croproi)
	for(i=1;i<m;i+=1)
		adja[croproi[i]][croproi[i-1]]=cle+1
		adja[croproi[i-1]][croproi[i]]=cle+1
		parent[croproi[i]]=croproi[0]
	endfor
endfor
//mettre à jour les interfaces
Duplicate/O calcDegree(adja) Degree_Obs
Duplicate/O calcChef(parent) Chef_Obs
Duplicate/O CompteTribue(chef_Obs,adja,wave0stripped,Parent) Tribue_Obs
Duplicate/O Sort2D(Tribue_Obs,1,1,-1) Tribue_Obs
Make/O/D/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphSel=0
Make/O/T/N=(dimsize(Chef_Obs,0),dimsize(ResGraphTitle,0)) ResGraphList=""
ResGraphList=num2str(Tribue_Obs)
Make/O/D/N=(n,3)/FREE Obs_Vizu
Obs_Vizu[][0]=wave0stripped[x]
Obs_Vizu[][1]=wave1stripped[x]
Obs_Vizu[][2]=wave2stripped[x]
Duplicate/O buildTrace(Adja,Obs_vizu,Parent), Traceur_obs
end



Function SpecluSet(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			speclu4chromato(dval)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc_7(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			Kmean4chromato(dval)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc_18(ba) : ButtonControl
	STRUCT WMButtonAction &ba
wave adja_obs,parent_obs,roipnts,Clu_Cla3
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			ReformateAdja4chromato(adja_obs,parent_obs,roipnts,Clu_Cla3,3)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function plotAllChromatograms()
wave/T advchromalist
variable n=dimsize(advchromalist,0),k
display/K=1
string lex,ley
for(k=0;k<n;k+=1)
	ley=advchromalist[k][1]+"_2chroma"
	lex=advchromalist[k][1]+"_3chroma"
	Duplicate/O/FREE $ley ord
	Duplicate/O/FREE $lex abc
	appendtograph $ley vs $lex
	ModifyGraph muloffset($ley)={0,1/wavemax(ord)}
	ModifyGraph offset($ley)={0,k}
endfor
end

function hayaFullReform(minSize)
variable minSize
wave/T ResGraphList
wave Parent_Obs,adja_obs,roipnts,Clu_Cla3
variable n=dimsize(ResGraphList,0),k
make/O/D/N=(n) chefs,taille
chefs=str2num(ResGraphList[x][0])
taille=str2num(ResGraphList[x][1])
for(k=0;k<n;k+=1)
	if(taille[k]>=minSize)
		Arbre2ROI(Parent_Obs,chefs[k])
		speClu4chromato(0)
		ReformateAdja4chromato(adja_obs,parent_obs,roipnts,Clu_Cla3,3)
		doupdate
	else
		k=n
	endif
endfor
end

function Qspan()
wave parent_obs
wave ciblebag_MurchMill
variable k, n=dimsize(ciblebag_MurchMill,0)
for(k=0;k<n;k+=1)
	Arbre2ROI(Parent_Obs,ciblebag_MurchMill[k][0])
	doupdate
endfor
end

function/WAVE MolbagPairTest(name1,name2)
string name1,name2

string sparent1="Parent_"+name1
string sparent2="Parent_"+name2
wave parent1=$sparent1
wave parent2=$sparent2

string windowname
string squada,squadb

string smolbag1="Molbag_"+name1
string schargebag1="Chargebag_"+name1
string slistbag1="Listbag_"+name1
string sciblebag1="Ciblebag_"+name1
wave molbag1=$smolbag1
wave chargebag1=$schargebag1
wave/T listbag1=$slistbag1
wave ciblebag1=$sciblebag1
string smolbag2="Molbag_"+name2
string schargebag2="Chargebag_"+name2
string slistbag2="Listbag_"+name2
string sciblebag2="Ciblebag_"+name2
wave molbag2=$smolbag2
wave chargebag2=$schargebag2
wave/T listbag2=$slistbag2
wave ciblebag2=$sciblebag2
variable n=dimsize(ciblebag1,0), m=dimsize(ciblebag2,0),k,j,c=0,tot=n*m
make/O/D/N=(tot,5)/FREE testRes
for(k=0;k<n;k+=1)
	for(j=0;j<m;j+=1)
		//N=1 et O=0 dans les deux
		testRes[c][0]=molbag1[6][k]==1 && molbag2[6][j]==1 && molbag1[7][k]==0 && molbag2[7][j]==0
		//DBE egal
		testRes[c][1]=(molbag1[5][k]+molbag1[113][k]+molbag1[112][k]+molbag1[6][k]/2-molbag1[0][k]/2)==(molbag2[5][j]+molbag2[113][j]+molbag2[112][j]+molbag2[6][j]/2-molbag2[0][j]/2)
		testRes[c][2]=(molbag1[5][k]+molbag1[113][k]+molbag1[112][k]+molbag1[6][k]/2-molbag1[0][k]/2)
		testRes[c][3]=ciblebag1[k][0]
		testRes[c][4]=ciblebag2[j][0]
		if(testRes[c][0] && testRes[c][1])
			windowname="DBE_equals_"+num2str(0.5+testRes[c][2])
			squada="Quad_"+name1+"_"+num2str(ciblebag1[k][0])
			squadb="Quad_"+name2+"_"+num2str(ciblebag2[j][0])
			if(wintype(windowname)==0)
				display/K=1/N=$windowname
			endif
			if(waveExists($squada)==0)
				Duplicate/O quadFromTree(parent1,ciblebag1[k][0],name1), $squada
				wave quada=$squada
				appendtograph/W=$windowname quada[][1] vs quada[][0]
				ModifyGraph mode($squada)=8,marker($squada)=19,useMrkStrokeRGB($squada)=1
				ModifyGraph zColor($squada)={$squada[][3],*,1500,ColdWarm,0}
				ModifyGraph muloffset($squada)={0,maxIntQuad($squada)}
				print squada,maxIntQuad($squada)
			endif
			if(waveExists($squadb)==0)
				Duplicate/O quadFromTree(parent2,ciblebag2[j][0],name2), $squadb
				wave quadb=$squadb
				appendtograph/W=$windowname quadb[][1] vs quadb[][0]
				ModifyGraph mode($squadb)=8,marker($squadb)=16,useMrkStrokeRGB($squadb)=1
				ModifyGraph zColor($squadb)={$squadb[][3],*,1500,ColdWarm,0}
				ModifyGraph muloffset($squadb)={0,maxIntQuad($squadb)}
				ModifyGraph offset($squadb)={0,1}
			endif
			//print c, testRes[c][2], testRes[c][3],testRes[c][4]
		endif
		c+=1
	endfor
endfor
return testRes
end

function killQuads()
string lesquads=wavelist("quad_*",";","")
string lequad=""
variable k,n=itemsInList(lesquads)
for(k=0;k<n;k+=1)
	lequad=stringfromList(k,lesquads)
	killwaves $lequad
endfor
end

function maxIntQuad(quad)
wave quad
extract/FREE quad,crop,y==1
return 1/wavemax(crop)
end

function modifyDBEequals()
string lesDBE=winlist("DBE_equals_*",";","")
string leDBE=""
string formula
variable k,n=itemsinlist(lesDBE),DBE
for(k=0;k<n;k+=1)
	leDBE=stringfromList(k,lesDBE)
	DBE=str2num(leDBE[11])
	SetAxis/W=$leDBE bottom 100,500;DelayUpdate
	ModifyGraph/W=$leDBE grid(left)=1,nticks(left)=2,fStyle=1,noLabel(left)=1,axThick=2,gridRGB(left)=(0,0,0);DelayUpdate
	ModifyGraph/W=$leDBE grid=1,nticks(bottom)=10,gridStyle(bottom)=4,gridRGB=(52428,52428,52428)
	Label/W=$leDBE left "Normalized intensity\r\\W5019Murchison\t\t\\W5016Ryugu";DelayUpdate
	Label/W=$leDBE bottom "Mass/Charge (u.)"
	formula="\Z18\f01C\Bn\M\Z18\f01H\B2n"+num2str(4-2*DBE)+"\M\Z18\f01N\S+"
	TextBox/W=$leDBE/K/N=text0
	TextBox/W=$leDBE/C/N=text0/A=MC formula
	TextBox/W=$leDBE/C/N=text0/A=RT/X=0/Y=0
	ModifyGraph/W=$leDBE axOffset(bottom)=1,lblMargin(bottom)=10
	movewindow/W=$leDBE 0,0,400,300
endfor
end

function/WAVE reshapeHeatmap(name)
string name
string Smass=name+"_4chroma"
string Sdata=name+"_5chroma"
wave mass=$Smass
wave data=$Sdata
variable n,m,k,j
Duplicate/O/FREE data res
res=0
n=dimsize(res,0)
m=dimsize(res,1)
setscale x,wavemin(mass),wavemax(mass), res
for(k=0;k<n;k+=1)
	for(j=0;j<m;j+=1)
		res[scaletoindex(res,mass[k],0)][j]+=data[k][j]
	endfor
endfor
return res
end

function externalization()
	print "oui ça marche"
end