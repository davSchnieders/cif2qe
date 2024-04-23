#!/usr/bin/python3
from ase import io
from ase import build
from ase.io.espresso import write_espresso_in
import sys
import os
import shutil
import click


def read_cif(filename):
    return io.read(filename)


def prepare_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def get_elements(atoms):
    seen = set()
    elements = [x for x in atoms.get_chemical_symbols() if x not in seen and not seen.add(x)]
    element_list = []
    for element in elements:
        element_list.append(
            {
                'symbol': element,
                'n_of_element': atoms.get_chemical_symbols().count(element)
            }
        )

    return element_list


def get_pseudo_potentials(atoms, max_valence, use_polarization):
    pseudo_potentials = {}
    pp_dir = os.environ['PP_DIR']
    elements = get_elements(atoms)
    for element in elements:
        files = []
        for file in os.listdir(pp_dir):
            if file.startswith('{element}.'.format(element=element['symbol'])):
                shells = []
                n_aos = 0
                n_el = 0
                with open("{pp_dir}/{file}".format(pp_dir=pp_dir, file=file), 'r') as pp_file:
                    while (True):
                        if pp_file.readline().strip() == "Valence configuration:":
                            break
                    pp_file.readline()
                    while (True):
                        line_raw = pp_file.readline()
                        if line_raw.strip() == "Generation configuration:":
                            break
                        line = line_raw.split()
                        if float(line[3]) > 0.0 or use_polarization:
                            shells.append(line[0])
                            n_aos += int(line[2]) * 2 + 1
                            n_el += int(float(line[3]))
                files.append({
                    'name': file,
                    'shells': ','.join(shells),
                    'n_aos': n_aos,
                    'n_el': n_el
                })

        i_pp_max_aos = -1
        max_aos = -1

        print("\n\nFound the following pp files for {element}".format(element=element))
        for i, file in enumerate(files):
            if file['n_aos'] > max_aos:
                i_pp_max_aos = i
                max_aos = file['n_aos']

            print(
                '{num}. {name} (val.: {shells}; n_aos: {n_aos}; n_el: {n_el})'.format(num=i + 1, name=file['name'],
                                                                       shells=file['shells'],
                                                                       n_aos=file['n_aos'],
                                                                       n_el=n_el))
        print('Which PP should we use? (Max AOs: {i_pp_max_aos})'.format(i_pp_max_aos=i_pp_max_aos+1))

        def get_input_check_within_range(n):
            try:
                inp = int(input())
            except ValueError:
                print("Invalid argument. Please provide a number up to {n}!".format(n=n))
                return get_input_check_within_range(n)
            if inp > n:
                print("Invalid argument. Please provide a number up to {n}!".format(n=n))
                return get_input_check_within_range(n)
            return inp - 1

        i_pp = i_pp_max_aos if max_valence else get_input_check_within_range(len(files))
        print('\nChoose {i_pp}!\n'.format(i_pp=i_pp+1))
        pp_to_use = files[i_pp]['name']

        shutil.copyfile('{pp_dir}/{file}'.format(pp_dir=pp_dir,file=pp_to_use),
                        '{dir}/{file}'.format(dir=str(atoms.symbols), file=pp_to_use))

        pseudo_potentials[element['symbol']] = {
            'file': pp_to_use,
            'n_aos': files[i_pp]['n_aos'],
            'n_el': files[i_pp]['n_el'],
            'shells': files[i_pp]['shells']
        }

    return pseudo_potentials


@click.command()
@click.option('--sc', default="1 1 1", help='Supercell Dimensions.')
@click.option('--metal', is_flag=True, default=False, help='Set if a metal is to be calculated')
@click.option('--max_valence', is_flag=True, default=False, help='If set, maximum valence will be used troughout.')
@click.option('--use_polarization', is_flag=True, default=False, help='If set, polarization functions will be used troughout.')
@click.argument('cif')
def run(cif, sc, metal, max_valence, use_polarization):
    sc_array = [int(i) for i in sc.split()]
    print(sc_array)
    atoms_raw = read_cif(cif)
    # if atoms_raw.get_cell().volume > 200:
    #     print("Cell too big!: {vol}".format(vol=atoms_raw.get_cell().volume))
    #     return 1
    atoms = build.make_supercell(atoms_raw, [[sc_array[0], 0, 0], [0, sc_array[1], 0], [0, 0, sc_array[2]]])
    elements = get_elements(atoms)
    prepare_dir(str(atoms.symbols))
    prepare_dir("{dir}/remove_sym".format(dir=str(atoms.symbols)))
    pseudo_potentials = get_pseudo_potentials(atoms, max_valence, use_polarization)

    shutil.copyfile(cif, "{dir}/{dir}.cif".format(dir=str(atoms.symbols)))

    pseudopotentials = {element: pp['file'] for element, pp in pseudo_potentials.items()}
    n_bands = sum([element['n_of_element'] * pseudo_potentials[element['symbol']]['n_aos'] for element in elements])
    n_el_tot = sum([element['n_of_element'] * pseudo_potentials[element['symbol']]['n_el'] for element in elements])

    input_data = {
        'control': {
            'calculation': 'scf',
            'prefix': str(atoms.symbols),
            'outdir': './',
            'wf_collect': True,
            'pseudo_dir': './'
        },
        'system': {
            'nosym': False,
            'ecutwfc': 60,
            'nbnd': n_bands,
            'nspin': 1 if n_el_tot % 2 == 0 else 2
        },
        'electrons': {
            'conv_thr': 1.0e-10,
            'electron_maxstep': 100
        }
    }

    if metal:
        input_data['system']['degauss'] = 0.01
        input_data['system']['occupations'] = 'smearing'
        input_data['system']['smearing'] = 'cold'

    with open('{dir}/pw.scf.in'.format(dir=str(atoms.symbols)), 'w') as pwfile:
        write_espresso_in(fd=pwfile, atoms=atoms, kspacing=0.025, input_data=input_data,
                          pseudopotentials=pseudopotentials,
                          crystal_coordinates=True)

    input_data['control']['calculation'] = 'nscf'
    input_data['system']['nosym'] = True
    input_data['electrons'] = {
        'startingpot': 'file'
    }

    with open('{dir}/remove_sym/pw.scf.in'.format(dir=str(atoms.symbols)), 'w') as pwfile:
        write_espresso_in(fd=pwfile, atoms=atoms, kspacing=0.025, input_data=input_data,
                          pseudopotentials=pseudopotentials,
                          crystal_coordinates=True)

    with open('{dir}/remove_sym/lobsterin'.format(dir=str(atoms.symbols)), 'w') as lobsterinfile:
        lobsterin_text = """
    cohpstartenergy -1
    cohpendenergy 1
    cohpsteps 2

    saveprojectiontofile
    skipmadelungenergy

    basisset pbevaspfit2015

    """

        for element, pp in pseudo_potentials.items():
            lobsterin_text += "basisfunctions {element} {aos}\n".format(element=element,
                                                                        aos=pp['shells'].replace(',', ' '))

        lobsterinfile.write(lobsterin_text)


if __name__ == "__main__":
    run()

