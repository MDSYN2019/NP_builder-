import numpy as np
import scipy
import matplotlib.pyplot as plt

#############################################################################
### BUILDING LAMMPS DATA FILE OF INDIVIDUAL ALKANETHIOL-GOLD NANOPARTICLE ###
### United Atom Model                                                     ###
### Written by Wenbin Li, M.I.T., November, 2013                          ###
### Modified by S. Y. Noh, 2019                                           ###
###                                                                       ###
#############################################################################
    
###  The program outputs LAMMPS data file "data.np" and CFG configuration file
###  "gold_np_thiol.cfg", which can be visualized by Ju Li's Atomeye program.
###  Atomeye is a free, fast visualization software. It can be downloaded from
###  http://li.mit.edu/Archive/Graphics/A/
    
a = 4.0782 # gold lattice constant
BoxSize = 100 # supercell edge length
H = np.concatenate([[BoxSize,0,0],[0,BoxSize,0],[0,0,BoxSize]]) # supercell matrix
natom_thiol = 9 # so total number of atoms per thiol will be 9
R=input_('radius of the sphere to put thiols (say = 30):')
#

# UAModelSingleNP.m:19
    ###  number of atoms per thiol in united atom model.
    ###  for octanethiol SC8H17, it will be coarse-grained as SC8
    natom_thiol=9
# UAModelSingleNP.m:24
    ###  put the thiol molecules initially on a sphere centered
###  on the gold nanocrystal. The sphere has radius R
  UAModelSingleNP.m:28
    ###  the azimuthal angle and polar angle are divided evenly
###  forming a grid, the thiol molecules are then put on these grids
###  total number of thiol molecules will be the product
###  of nphi and ntheta
    nphi=8
# UAModelSingleNP.m:34
    ntheta=17
# UAModelSingleNP.m:35
    #### total number of thiols
    n_thiol=dot(nphi,ntheta)
# UAModelSingleNP.m:38
    ##########################################
### BUILDING ICOSAHEDRAL GOLD NANOCRYSTAL
###########################################
    
    ###  layer of gold atom in gold icosahedral nanocrystal
###  different number of gold layers will result in
###  different number of gold atoms in the nanocrystal
###  gold_layer = 6 results in a nanoparticle with 561 gold
###  atoms. The nanoparticle diameter will be around 3 nm
    
    gold_layer=6
# UAModelSingleNP.m:50
    n_gold=(dot(10,(gold_layer ** 3)) + dot(11,gold_layer)) / 3 - dot(5,(gold_layer ** 2)) - 1
# UAModelSingleNP.m:51
    ###  matrix store gold atom coordinates
    Pos_gold=zeros(n_gold,3)
# UAModelSingleNP.m:54
    ###  center of gold nanoparticle
    center=concat([0.5,0.5,0.5])
# UAModelSingleNP.m:57
    ###  parater used to build icosohedron
    phi=(sqrt(5) + 1) / 2
# UAModelSingleNP.m:60
    ### coordinates of the twelve vertices of icosahedron
    V=zeros(12,3)
# UAModelSingleNP.m:63
    V[1,arange()]=concat([1,0,phi])
# UAModelSingleNP.m:65
    V[2,arange()]=concat([- 1,0,phi])
# UAModelSingleNP.m:65
    V[3,arange()]=concat([0,- phi,1])
# UAModelSingleNP.m:65
    V[4,arange()]=concat([phi,- 1,0])
# UAModelSingleNP.m:66
    V[5,arange()]=concat([phi,1,0])
# UAModelSingleNP.m:66
    V[6,arange()]=concat([0,phi,1])
# UAModelSingleNP.m:66
    for i in arange(7,12).reshape(-1):
        V[i,arange()]=- V(i - 6,arange())
# UAModelSingleNP.m:69
    
    ### end points coordinates of the thirty eges, each edge has two end points
    E=zeros(2,3,30)
# UAModelSingleNP.m:73
    E[1,arange(),1]=concat([1,0,phi])
# UAModelSingleNP.m:75
    E[2,arange(),1]=concat([- 1,0,phi])
# UAModelSingleNP.m:75
    E[1,arange(),2]=concat([1,0,phi])
# UAModelSingleNP.m:76
    E[2,arange(),2]=concat([0,- phi,1])
# UAModelSingleNP.m:76
    E[1,arange(),3]=concat([1,0,phi])
# UAModelSingleNP.m:77
    E[2,arange(),3]=concat([phi,- 1,0])
# UAModelSingleNP.m:77
    E[1,arange(),4]=concat([1,0,phi])
# UAModelSingleNP.m:78
    E[2,arange(),4]=concat([phi,1,0])
# UAModelSingleNP.m:78
    E[1,arange(),5]=concat([1,0,phi])
# UAModelSingleNP.m:79
    E[2,arange(),5]=concat([0,phi,1])
# UAModelSingleNP.m:79
    E[1,arange(),6]=concat([0,phi,1])
# UAModelSingleNP.m:80
    E[2,arange(),6]=concat([phi,1,0])
# UAModelSingleNP.m:80
    E[1,arange(),7]=concat([0,phi,1])
# UAModelSingleNP.m:81
    E[2,arange(),7]=concat([0,phi,- 1])
# UAModelSingleNP.m:81
    E[1,arange(),8]=concat([0,phi,1])
# UAModelSingleNP.m:82
    E[2,arange(),8]=concat([- phi,1,0])
# UAModelSingleNP.m:82
    E[1,arange(),9]=concat([0,phi,1])
# UAModelSingleNP.m:83
    E[2,arange(),9]=concat([- 1,0,phi])
# UAModelSingleNP.m:83
    E[1,arange(),10]=concat([phi,1,0])
# UAModelSingleNP.m:84
    E[2,arange(),10]=concat([phi,- 1,0])
# UAModelSingleNP.m:84
    E[1,arange(),11]=concat([phi,1,0])
# UAModelSingleNP.m:85
    E[2,arange(),11]=concat([1,0,- phi])
# UAModelSingleNP.m:85
    E[1,arange(),12]=concat([phi,1,0])
# UAModelSingleNP.m:86
    E[2,arange(),12]=concat([0,phi,- 1])
# UAModelSingleNP.m:86
    E[1,arange(),13]=concat([phi,- 1,0])
# UAModelSingleNP.m:87
    E[2,arange(),13]=concat([1,0,- phi])
# UAModelSingleNP.m:87
    E[1,arange(),14]=concat([1,0,- phi])
# UAModelSingleNP.m:88
    E[2,arange(),14]=concat([0,phi,- 1])
# UAModelSingleNP.m:88
    E[1,arange(),15]=concat([0,phi,- 1])
# UAModelSingleNP.m:89
    E[2,arange(),15]=concat([- phi,1,0])
# UAModelSingleNP.m:89
    for i in arange(16,30).reshape(-1):
        E[arange(),arange(),i]=- E(arange(),arange(),i - 15)
# UAModelSingleNP.m:92
    
    ### coordinates of the vertices of twenty faces
    
    F=zeros(3,3,20)
# UAModelSingleNP.m:97
    F[1,arange(),1]=concat([1,0,phi])
# UAModelSingleNP.m:99
    F[2,arange(),1]=concat([- 1,0,phi])
# UAModelSingleNP.m:99
    F[3,arange(),1]=concat([0,- phi,1])
# UAModelSingleNP.m:99
    F[1,arange(),2]=concat([1,0,phi])
# UAModelSingleNP.m:100
    F[2,arange(),2]=concat([0,- phi,1])
# UAModelSingleNP.m:100
    F[3,arange(),2]=concat([phi,- 1,0])
# UAModelSingleNP.m:100
    F[1,arange(),3]=concat([1,0,phi])
# UAModelSingleNP.m:101
    F[2,arange(),3]=concat([phi,- 1,0])
# UAModelSingleNP.m:101
    F[3,arange(),3]=concat([phi,1,0])
# UAModelSingleNP.m:101
    F[1,arange(),4]=concat([1,0,phi])
# UAModelSingleNP.m:102
    F[2,arange(),4]=concat([phi,1,0])
# UAModelSingleNP.m:102
    F[3,arange(),4]=concat([0,phi,1])
# UAModelSingleNP.m:102
    F[1,arange(),5]=concat([1,0,phi])
# UAModelSingleNP.m:103
    F[2,arange(),5]=concat([0,phi,1])
# UAModelSingleNP.m:103
    F[3,arange(),5]=concat([- 1,0,phi])
# UAModelSingleNP.m:103
    F[1,arange(),6]=concat([0,phi,1])
# UAModelSingleNP.m:104
    F[2,arange(),6]=concat([phi,1,0])
# UAModelSingleNP.m:104
    F[3,arange(),6]=concat([0,phi,- 1])
# UAModelSingleNP.m:104
    F[1,arange(),7]=concat([0,phi,1])
# UAModelSingleNP.m:105
    F[2,arange(),7]=concat([0,phi,- 1])
# UAModelSingleNP.m:105
    F[3,arange(),7]=concat([- phi,1,0])
# UAModelSingleNP.m:105
    F[1,arange(),8]=concat([0,phi,1])
# UAModelSingleNP.m:106
    F[2,arange(),8]=concat([- phi,1,0])
# UAModelSingleNP.m:106
    F[3,arange(),8]=concat([- 1,0,phi])
# UAModelSingleNP.m:106
    F[1,arange(),9]=concat([phi,1,0])
# UAModelSingleNP.m:107
    F[2,arange(),9]=concat([phi,- 1,0])
# UAModelSingleNP.m:107
    F[3,arange(),9]=concat([1,0,- phi])
# UAModelSingleNP.m:107
    F[1,arange(),10]=concat([phi,1,0])
# UAModelSingleNP.m:108
    F[2,arange(),10]=concat([1,0,- phi])
# UAModelSingleNP.m:108
    F[3,arange(),10]=concat([0,phi,- 1])
# UAModelSingleNP.m:108
    for i in arange(11,20).reshape(-1):
        F[arange(),arange(),i]=- F(arange(),arange(),i - 10)
# UAModelSingleNP.m:111
    
    ### absolute coordinates of the inner-most gold layer
    scaling=a / sqrt(dot(2,(1 + dot(phi,phi))))
# UAModelSingleNP.m:115
    V=dot(V,scaling)
# UAModelSingleNP.m:117
    E=dot(E,scaling)
# UAModelSingleNP.m:118
    F=dot(F,scaling)
# UAModelSingleNP.m:119
    ### coordinates of the centeral atom
    Pos_gold[1,arange()]=center
# UAModelSingleNP.m:122
    ### coordinates of the second innermost layer atoms
    
    for i in arange(1,12).reshape(-1):
        Pos_gold[i + 1,1]=V(i,1) / BoxSize + center(1)
# UAModelSingleNP.m:127
        Pos_gold[i + 1,2]=V(i,2) / BoxSize + center(2)
# UAModelSingleNP.m:128
        Pos_gold[i + 1,3]=V(i,3) / BoxSize + center(3)
# UAModelSingleNP.m:129
    
    ### coordinates of the third innermost layer atoms
    for i in arange(1,12).reshape(-1):
        Pos_gold[i + 13,1]=dot(V(i,1),2) / BoxSize + center(1)
# UAModelSingleNP.m:134
        Pos_gold[i + 13,2]=dot(V(i,2),2) / BoxSize + center(2)
# UAModelSingleNP.m:135
        Pos_gold[i + 13,3]=dot(V(i,3),2) / BoxSize + center(3)
# UAModelSingleNP.m:136
    
    ### coordinates of the fourth innermost layer atoms
    for i in arange(1,30).reshape(-1):
        Pos_gold[i + 25,1]=(E(1,1,i) + E(2,1,i)) / BoxSize + center(1)
# UAModelSingleNP.m:141
        Pos_gold[i + 25,2]=(E(1,2,i) + E(2,2,i)) / BoxSize + center(2)
# UAModelSingleNP.m:142
        Pos_gold[i + 25,3]=(E(1,3,i) + E(2,3,i)) / BoxSize + center(3)
# UAModelSingleNP.m:143
    
    ### the other layers of atoms, which have vertice atoms, edges atoms and face atoms
    
    for N in arange(4,gold_layer).reshape(-1):
        #total atom up to N-th layer
        NatomToNow=(dot(10,((N - 1) ** 3)) + dot(11,(N - 1))) / 3 - dot(5,((N - 1) ** 2)) - 1
# UAModelSingleNP.m:150
        for i in arange(1,12).reshape(-1):
            Pos_gold[NatomToNow + i,1]=dot(V(i,1),(N - 1)) / BoxSize + center(1)
# UAModelSingleNP.m:154
            Pos_gold[NatomToNow + i,2]=dot(V(i,2),(N - 1)) / BoxSize + center(2)
# UAModelSingleNP.m:155
            Pos_gold[NatomToNow + i,3]=dot(V(i,3),(N - 1)) / BoxSize + center(3)
# UAModelSingleNP.m:156
        NatomToNow=NatomToNow + 12
# UAModelSingleNP.m:159
        for i in arange(1,30).reshape(-1):
            for j in arange(1,(N - 2)).reshape(-1):
                Pos_gold[NatomToNow + dot((i - 1),(N - 2)) + j,1]=(dot(E(1,1,i),(N - 1)) + dot((E(2,1,i) - E(1,1,i)),j)) / BoxSize + center(1)
# UAModelSingleNP.m:164
                Pos_gold[NatomToNow + dot((i - 1),(N - 2)) + j,2]=(dot(E(1,2,i),(N - 1)) + dot((E(2,2,i) - E(1,2,i)),j)) / BoxSize + center(2)
# UAModelSingleNP.m:165
                Pos_gold[NatomToNow + dot((i - 1),(N - 2)) + j,3]=(dot(E(1,3,i),(N - 1)) + dot((E(2,3,i) - E(1,3,i)),j)) / BoxSize + center(3)
# UAModelSingleNP.m:166
        NatomToNow=NatomToNow + dot(30,(N - 2))
# UAModelSingleNP.m:170
        for i in arange(1,20).reshape(-1):
            for m in arange(1,(N - 3)).reshape(-1):
                for n in arange(1,(N - m - 2)).reshape(-1):
                    Pos_gold[NatomToNow + dot((dot(2,N) - m - 4),(m - 1)) / 2 + n,1]=(dot(F(1,1,i),(N - 1)) + dot((dot(F(2,1,i),(N - 1)) - dot(F(1,1,i),(N - 1))),m) / (N - 1) + dot((dot(F(3,1,i),(N - 1)) - dot(F(1,1,i),(N - 1))),n) / (N - 1)) / BoxSize + center(1)
# UAModelSingleNP.m:176
                    Pos_gold[NatomToNow + dot((dot(2,N) - m - 4),(m - 1)) / 2 + n,2]=(dot(F(1,2,i),(N - 1)) + dot((dot(F(2,2,i),(N - 1)) - dot(F(1,2,i),(N - 1))),m) / (N - 1) + dot((dot(F(3,2,i),(N - 1)) - dot(F(1,2,i),(N - 1))),n) / (N - 1)) / BoxSize + center(2)
# UAModelSingleNP.m:178
                    Pos_gold[NatomToNow + dot((dot(2,N) - m - 4),(m - 1)) / 2 + n,3]=(dot(F(1,3,i),(N - 1)) + dot((dot(F(2,3,i),(N - 1)) - dot(F(1,3,i),(N - 1))),m) / (N - 1) + dot((dot(F(3,3,i),(N - 1)) - dot(F(1,3,i),(N - 1))),n) / (N - 1)) / BoxSize + center(3)
# UAModelSingleNP.m:180
            NatomToNow=NatomToNow + dot((N - 2),(N - 3)) / 2
# UAModelSingleNP.m:184
    
    ###################################################
### BUIDING ICOSAHEDRON GOLD NANOPARTICLES FINISHED
### NEXT BUILD THIOL MOLECULES
###################################################
    
    ### coordinates of the thiol molecules
    Pos_chain=zeros(natom_thiol,3,n_thiol)
# UAModelSingleNP.m:194
    ### matrix for conversion from absolute to reduced coordinates
    H1=concat([[1 / BoxSize,0,0],[0,1 / BoxSize,0],[0,0,1 / BoxSize]])
# UAModelSingleNP.m:197
    ###  the following vectors are used for determing the coordinates of thiol molecules
    ShiftVector1=concat([1.005,0,1.51])
# UAModelSingleNP.m:200
    ShiftVector2=concat([- 0.919,0,1.224])
# UAModelSingleNP.m:201
    TranslateVector=concat([- 0.071,0,2.49])
# UAModelSingleNP.m:202
    n_unit=natom_thiol - 1
# UAModelSingleNP.m:204
    k=1
# UAModelSingleNP.m:206
    ###  increment of azimuthal and polar angles
    dphi=dot(2,pi) / nphi
# UAModelSingleNP.m:209
    dtheta=pi / (ntheta + 3)
# UAModelSingleNP.m:209
    for phi in arange(dphi,dot(2,pi),dphi).reshape(-1):
        for theta in arange(dot(2,dtheta),dot(dtheta,(ntheta + 1)),dtheta).reshape(-1):
            Pos_chain[1,1,k]=BoxSize / 2 + dot(dot(R,sin(theta)),cos(phi))
# UAModelSingleNP.m:213
            Pos_chain[1,2,k]=BoxSize / 2 + dot(dot(R,sin(theta)),sin(phi))
# UAModelSingleNP.m:214
            Pos_chain[1,3,k]=BoxSize / 2 + dot(R,cos(theta))
# UAModelSingleNP.m:215
            M=concat([dot(cos(theta),cos(phi)),dot(cos(theta),sin(phi)),- sin(theta),- sin(phi),cos(phi),0,dot(sin(theta),cos(phi)),dot(sin(theta),sin(phi)),cos(theta)])
# UAModelSingleNP.m:218
            Pos_chain[2,arange(),k]=Pos_chain(1,arange(),k) + dot(ShiftVector1,M)
# UAModelSingleNP.m:223
            Pos_chain[3,arange(),k]=Pos_chain(1,arange(),k) + dot((ShiftVector1 + ShiftVector2),M)
# UAModelSingleNP.m:224
            for n in arange(1,(n_unit - 2) / 2).reshape(-1):
                Pos_chain[3 + dot(2,n) - 1,arange(),k]=Pos_chain(1,arange(),k) + dot((ShiftVector1 + dot(TranslateVector,n)),M)
# UAModelSingleNP.m:227
                Pos_chain[3 + dot(2,n),arange(),k]=Pos_chain(1,arange(),k) + dot((ShiftVector1 + ShiftVector2 + dot(TranslateVector,n)),M)
# UAModelSingleNP.m:228
            ### Conversion from absolute coordinates to reduced coordinates
            Pos_chain[arange(),arange(),k]=dot(Pos_chain(arange(),arange(),k),H1)
# UAModelSingleNP.m:232
            k=k + 1
# UAModelSingleNP.m:234
    
    ##########################################################################
    
    ### COORDINATES DETERMINATION ENDS HERE, NEXT OUTPUT CFG COORDINATION FILE
### THE CFG COORDINATION FILE CAN BE VISUALIZED BY JU LI'S ATOMEYE PROGRAM
### http://li.mit.edu/Archive/Graphics/A/
    
    ##########################################################################
    
    atom_name_Au='Au'
# UAModelSingleNP.m:246
    atom_mass_Au=196.96655
# UAModelSingleNP.m:247
    atom_name_S='S'
# UAModelSingleNP.m:249
    atom_mass_S=32.065
# UAModelSingleNP.m:250
    atom_name_C='C'
# UAModelSingleNP.m:252
    atom_mass_C=12.001
# UAModelSingleNP.m:253
    natom_total=n_gold + dot(natom_thiol,n_thiol)
# UAModelSingleNP.m:255
    cfg_name='gold_np_thiol.cfg'
# UAModelSingleNP.m:257
    cfg=fopen(cfg_name,'w')
# UAModelSingleNP.m:258
    fprintf(cfg,'Number of particles = %d\n',natom_total)
    fprintf(cfg,'H0(1,1)= %.30g A\n',H(1,1))
    fprintf(cfg,'H0(1,2)= %.30g A\n',H(1,2))
    fprintf(cfg,'H0(1,3)= %.30g A\n',H(1,3))
    fprintf(cfg,'H0(2,1)= %.30g A\n',H(2,1))
    fprintf(cfg,'H0(2,2)= %.30g A\n',H(2,2))
    fprintf(cfg,'H0(2,3)= %.30g A\n',H(2,3))
    fprintf(cfg,'H0(3,1)= %.30g A\n',H(3,1))
    fprintf(cfg,'H0(3,2)= %.30g A\n',H(3,2))
    fprintf(cfg,'H0(3,3)= %.30g A\n',H(3,3))
    fprintf(cfg,'.NO_VELOCITY.\n')
    fprintf(cfg,'entry_count = 3\n')
    fprintf(cfg,'%g\n',atom_mass_Au)
    fprintf(cfg,'%2s\n',atom_name_Au)
    for i in arange(1,n_gold).reshape(-1):
        fprintf(cfg,'%.30g %.30g %.30g\n',Pos_gold(i,1),Pos_gold(i,2),Pos_gold(i,3))
    
    fprintf(cfg,'%g\n',atom_mass_S)
    fprintf(cfg,'%2s\n',atom_name_S)
    for k in arange(1,n_thiol).reshape(-1):
        fprintf(cfg,'%.30g %.30g %.30g\n',Pos_chain(1,1,k),Pos_chain(1,2,k),Pos_chain(1,3,k))
    
    
    fprintf(cfg,'%g\n',atom_mass_C)
    fprintf(cfg,'%2s\n',atom_name_C)
    for k in arange(1,n_thiol).reshape(-1):
        for i in arange(2,natom_thiol).reshape(-1):
            fprintf(cfg,'%.30g %.30g %.30g\n',Pos_chain(i,1,k),Pos_chain(i,2,k),Pos_chain(i,3,k))
    
    
    fclose(cfg)
    #convert to absolute coordinates
    
    Pos_gold=dot(Pos_gold,H)
# UAModelSingleNP.m:300
    for k in arange(1,n_thiol).reshape(-1):
        Pos_chain[arange(),arange(),k]=dot(Pos_chain(arange(),arange(),k),H)
# UAModelSingleNP.m:303
    
    ##############################################################
    
    ###  NEXT OUTPUT LAMMPS DATA FILE
    
    ##############################################################
    
    ###  total number of atoms
    total_atoms=n_gold + dot(natom_thiol,n_thiol)
# UAModelSingleNP.m:315
    ###  total number of bonds
    total_bonds=dot((natom_thiol - 1),n_thiol)
# UAModelSingleNP.m:318
    ###  total number of angles
    total_angles=dot((natom_thiol - 2),n_thiol)
# UAModelSingleNP.m:321
    ###  total nubmer of dihedrals
    total_dihedrals=dot((natom_thiol - 3),n_thiol)
# UAModelSingleNP.m:324
    ###  atoms types:
###  Au-1
###  S-2
###  C-3, from coarse graining of CH2
###  C-4, from coarse graining of CH3
    atom_types=4
# UAModelSingleNP.m:331
    ###  bond type
###  1  S-C
###  2  C-C between CH2 and CH2
###  3  C-C between CH2 and CH3
    bond_types=3
# UAModelSingleNP.m:337
    ###  angle types
###  1  S-CH2-CH2
###  2  CH2-CH2-CH2
###  3  CH2-CH2-CH3
    angle_types=3
# UAModelSingleNP.m:343
    ###  dihedral types
###  1 S-CH2-CH2-CH2
###  2 CH2-CH2-CH2-CH2
###  3 CH2-CH2-CH2-CH3
    dihedral_types=3
# UAModelSingleNP.m:349
    atom_mass_Au=196.9666
# UAModelSingleNP.m:352
    atom_mass_S=32.065
# UAModelSingleNP.m:353
    atom_mass_CH2=14.0026
# UAModelSingleNP.m:354
    atom_mass_CH3=15.0345
# UAModelSingleNP.m:355
    ###  file name
    
    file_name='data.np'
# UAModelSingleNP.m:359
    fp=fopen(file_name,'w')
# UAModelSingleNP.m:361
    fprintf(fp,'Gold Nanocrystal - Alkanethiol System\n')
    fprintf(fp,'\n')
    fprintf(fp,'        %d   atoms\n',total_atoms)
    fprintf(fp,'        %d   bonds\n',total_bonds)
    fprintf(fp,'        %d   angles\n',total_angles)
    fprintf(fp,'        %d   dihedrals\n',total_dihedrals)
    fprintf(fp,'\n')
    fprintf(fp,'        %d   atom types\n',atom_types)
    fprintf(fp,'        %d   bond types\n',bond_types)
    fprintf(fp,'        %d   angle types\n',angle_types)
    fprintf(fp,'        %d   dihedral types\n',dihedral_types)
    fprintf(fp,'\n')
    fprintf(fp,'  0.00  %.5g        xlo xhi\n',H(1,1))
    fprintf(fp,'  0.00  %.5g        ylo yhi\n',H(2,2))
    fprintf(fp,'  0.00  %.5g        zlo zhi\n',H(3,3))
    fprintf(fp,'\n')
    fprintf(fp,'Masses\n')
    fprintf(fp,'\n')
    fprintf(fp,' 1 %g\n',atom_mass_Au)
    fprintf(fp,' 2 %g\n',atom_mass_S)
    fprintf(fp,' 3 %g\n',atom_mass_CH2)
    fprintf(fp,' 4 %g\n',atom_mass_CH3)
    #####################
###  PRINT ATOMS  ###
#####################
    
    fprintf(fp,'\n')
    fprintf(fp,'Atoms\n')
    fprintf(fp,'\n')
    ### PRINT FORMAT: ATOM_ID MOLECULE_ID ATOM_TYPE POSITION_X POSITON_Y POSITION_Z
    
    ### print gold atoms
    for i in arange(1,n_gold).reshape(-1):
        fprintf(fp,'   %d     1   1   %.5g  %.5g  %.5g\n',i,Pos_gold(i,1),Pos_gold(i,2),Pos_gold(i,3))
    
    for k in arange(1,n_thiol).reshape(-1):
        ### print sulfur atoms
        fprintf(fp,'   %d     %d   2  %.5g  %.5g  %.5g\n',n_gold + dot((k - 1),natom_thiol) + 1,1 + k,Pos_chain(1,1,k),Pos_chain(1,2,k),Pos_chain(1,3,k))
        for i in arange(2,(natom_thiol - 1)).reshape(-1):
            fprintf(fp,'   %d     %d   3  %.5g  %.5g  %.5g\n',n_gold + dot((k - 1),natom_thiol) + i,1 + k,Pos_chain(i,1,k),Pos_chain(i,2,k),Pos_chain(i,3,k))
        ### print coarse-grained CH3
        fprintf(fp,'   %d     %d   4  %.5g  %.5g  %.5g\n',n_gold + dot(k,natom_thiol),1 + k,Pos_chain(natom_thiol,1,k),Pos_chain(natom_thiol,2,k),Pos_chain(natom_thiol,3,k))
    
    #####################
###  PRINT BONDS  ###
#####################
    
    fprintf(fp,'\n')
    fprintf(fp,'Bonds\n')
    fprintf(fp,'\n')
    ### PRINT FORMAT: BOND_NUMBER BOND_TYPE ATOM_ID1 ATOM_ID2
    
    for k in arange(1,n_thiol).reshape(-1):
        ### print S-C bonds
        fprintf(fp,'     %d   1     %d     %d\n',dot((k - 1),(natom_thiol - 1)) + 1,n_gold + dot((k - 1),natom_thiol) + 1,n_gold + dot((k - 1),natom_thiol) + 2)
        for i in arange(2,(natom_thiol - 2)).reshape(-1):
            fprintf(fp,'     %d   2     %d     %d\n',dot((k - 1),(natom_thiol - 1)) + i,n_gold + dot((k - 1),natom_thiol) + i,n_gold + dot((k - 1),natom_thiol) + i + 1)
        ### print CH2-CH3 bonds
        fprintf(fp,'     %d   3     %d     %d\n',dot(k,(natom_thiol - 1)),n_gold + dot(k,natom_thiol) - 1,n_gold + dot(k,natom_thiol))
    
    ######################
###  Print Angles  ###
######################
    
    fprintf(fp,'\n')
    fprintf(fp,'Angles\n')
    fprintf(fp,'\n')
    ### PRINT FORMAT: ANGLE_NUMBER ANGLE_TYPE ATOM_ID1 ATOM_ID2 ATOM_ID3
    
    for k in arange(1,n_thiol).reshape(-1):
        ### S-CH2-CH2
        fprintf(fp,'     %d   1     %d     %d     %d\n',dot((k - 1),(natom_thiol - 2)) + 1,n_gold + dot((k - 1),natom_thiol) + 1,n_gold + dot((k - 1),natom_thiol) + 2,n_gold + dot((k - 1),natom_thiol) + 3)
        for i in arange(2,(natom_thiol - 3)).reshape(-1):
            fprintf(fp,'     %d   2     %d     %d     %d\n',dot((k - 1),(natom_thiol - 2)) + i,n_gold + dot((k - 1),natom_thiol) + i,n_gold + dot((k - 1),natom_thiol) + i + 1,n_gold + dot((k - 1),natom_thiol) + i + 2)
        ### CH2-CH2-CH3
        fprintf(fp,'     %d   3     %d     %d     %d\n',dot(k,(natom_thiol - 2)),n_gold + dot(k,natom_thiol) - 2,n_gold + dot(k,natom_thiol) - 1,n_gold + dot(k,natom_thiol))
    
    ##########################
###   PRINT DIHEDRALS  ###
##########################
    
    fprintf(fp,'\n')
    fprintf(fp,'Dihedrals\n')
    fprintf(fp,'\n')
    ### PRINT FORMAT: DIHEDRAL_N_GOLD DIHEDRAL_TYPE ATOM_ID1 ATOM_ID2 ATOM_ID3
    
    for k in arange(1,n_thiol).reshape(-1):
        ### S-CH2-CH2-CH2
        fprintf(fp,'     %d   1     %d     %d     %d     %d\n',dot((k - 1),(natom_thiol - 3)) + 1,n_gold + dot((k - 1),natom_thiol) + 1,n_gold + dot((k - 1),natom_thiol) + 2,n_gold + dot((k - 1),natom_thiol) + 3,n_gold + dot((k - 1),natom_thiol) + 4)
        for i in arange(2,(natom_thiol - 4)).reshape(-1):
            fprintf(fp,'     %d   2     %d     %d     %d     %d\n',dot((k - 1),(natom_thiol - 3)) + i,n_gold + dot((k - 1),natom_thiol) + i,n_gold + dot((k - 1),natom_thiol) + i + 1,n_gold + dot((k - 1),natom_thiol) + i + 2,n_gold + dot((k - 1),natom_thiol) + i + 3)
        ### CH2-CH2-CH2-CH3
        fprintf(fp,'     %d   3     %d     %d     %d     %d\n',dot(k,(natom_thiol - 3)),n_gold + dot(k,natom_thiol) - 3,n_gold + dot(k,natom_thiol) - 2,n_gold + dot(k,natom_thiol) - 1,n_gold + dot(k,natom_thiol))
    
    fclose(fp)
