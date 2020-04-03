from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
from plotting import energy_profile_plot
# -*- coding: utf-8 -*-

def optimization(pdb_file,options):
    """This function needs as an input a PDB file (the model), and gives as an output the optimized model."""
    env = environ()
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    path, ext = pdb_file.split('.')
    list = path.split('/')
    #if results/4g83/4g83model.pdb list[0]=results list[1]=4g83 list[2]=4g83model
    dir=list[0]+"/"+list[1]
    code=list[2]

    mdl = complete_pdb(env, pdb_file)
    mdl.write(file=path+'.ini')

    atmsel = selection(mdl)

    mpdf_ini = atmsel.energy()
    z_score_ini = mdl.assess_normalized_dope()
    mdl_ep_ini = atmsel.get_dope_profile()
    mdl_ep_ini_smoothed = mdl_ep_ini.get_smoothed()
    energy_profile_txt_path = dir + '/' + code + '_DOPE_EnergyProfile.txt'
    mdl_ep_ini_smoothed.write_to_file(energy_profile_txt_path)
    print("The unoptimized model's energy of " + code + " is: " + str(mpdf_ini[0]))
    print("The unoptimized Z-score of " + code + " is: " + str(z_score_ini))

    energy_profile_txt_path_opt=None

    if options.optimize:
        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trcfil = open(path+'.D00000001', 'w')
        cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))
        md.optimize(atmsel, temperature=300, max_iterations=50, actions=[actions.write_structure(10, path+'.D9999%04d.pdb'),actions.trace(10, trcfil)])
        cg.optimize(atmsel, max_iterations=20, actions=[actions.trace(5, trcfil)])
        mpdf = atmsel.energy()
        z_score = mdl.assess_normalized_dope()
        print("The final model energy of " + path + " is " + str(mpdf[0]))
        print("The final model energy of " + path + " is " + str(mpdf_ini[0]))
        print("The final z-score of " + code + " is: " + str(z_score))

        mdl.write(file=path+'_optimized.pdb')

        mdl_final = atmsel.get_dope_profile()
        mdl_final_smoothed = mdl_final.get_smoothed(window=50)
        energy_profile_txt_path_opt = dir + '/' + code + '_optimized_DOPE_EnergyProfile.txt'
        mdl_final_smoothed.write_to_file(energy_profile_txt_path_opt)
        mdl.write(file=dir + '/' + code + '_optimized.pdb')

    energy_profile_plot(options, dir, code, energy_profile_txt_path, energy_profile_txt_path_opt)
