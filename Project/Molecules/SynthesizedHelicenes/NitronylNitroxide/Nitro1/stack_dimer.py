from openbabel import openbabel

def stack_dimer(input_mol_file, stack_distance):

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mol")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_mol_file)

    mol_copy = openbabel.OBMol(mol)
    for atom in openbabel.OBMolAtomIter(mol_copy):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        atom.SetVector(x, y, z + stack_distance)

    mol += mol_copy
    mol.PerceiveBondOrders()  

    output_mol_file = "stacked_dimer_Nitro2.mol"
    obConversion.WriteFile(mol, output_mol_file)
    
    return output_mol_file

if __name__ == "__main__":
    input_mol_file = "opt_radical_Nitro2.mol"
    stacked_dimer_file = stack_dimer(input_mol_file, 5)
    print(f"Stacked dimer file created: {stacked_dimer_file}")
