"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    Phosphatidylglycerols consist of a glycerol backbone with fatty acyl chains at sn-1 and sn-2 positions,
    and a phosphatidyl group attached to another glycerol at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify phosphate groups in the molecule
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphate_atoms:
        return False, "No phosphate group found"

    # Iterate over phosphate groups to find potential phosphatidylglycerol structure
    for phospho_atom in phosphate_atoms:
        # Get oxygen atoms attached to phosphate
        phosphate_oxygens = [nbr for nbr in phospho_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(phosphate_oxygens) < 3:
            continue  # Not a full phosphate group
        
        # Initialize variables to track connected groups
        diacylglycerol_found = False
        headgroup_glycerol_found = False

        # Check oxygens connected to phosphate
        for oxy in phosphate_oxygens:
            # Get the atoms connected to the oxygen (excluding phosphate)
            oxy_neighbors = [nbr for nbr in oxy.GetNeighbors() if nbr.GetIdx() != phospho_atom.GetIdx()]
            if not oxy_neighbors:
                continue  # No other connections
            neighbor_atom = oxy_neighbors[0]
            
            if neighbor_atom.GetAtomicNum() == 6:
                # Check if this carbon is part of a glycerol backbone
                # Find three-carbon chain starting from this carbon
                paths = Chem.FindAllPathsOfLengthN(mol, 3, useBonds=True, useHs=False)
                for path in paths:
                    if neighbor_atom.GetIdx() in path:
                        # Get the atoms in the path
                        atoms_in_path = [mol.GetAtomWithIdx(idx) for idx in path]
                        # Check if all atoms are carbons
                        if all(atom.GetAtomicNum() == 6 for atom in atoms_in_path):
                            # Check if each carbon has an oxygen attached
                            oxygens_attached = []
                            for atom in atoms_in_path:
                                oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                                if not oxygens:
                                    break  # No oxygen attached
                                oxygens_attached.append(oxygens)
                            else:
                                # Check substitutions at sn-1, sn-2, sn-3 positions
                                # sn-1 and sn-2 positions should have ester bonds to fatty acids
                                sn1_carbon = atoms_in_path[0]
                                sn2_carbon = atoms_in_path[1]
                                sn3_carbon = atoms_in_path[2]
                                
                                # Check for ester linkage at sn-1 and sn-2
                                def is_ester_linked(carbon_atom):
                                    for oxygen_atom in [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]:
                                        for neighbor in oxygen_atom.GetNeighbors():
                                            if neighbor.GetIdx() == carbon_atom.GetIdx():
                                                continue
                                            if neighbor.GetAtomicNum() == 6:
                                                # Check if neighbor is carbonyl carbon (C=O)
                                                for bond in neighbor.GetBonds():
                                                    if bond.GetOtherAtomIdx(neighbor.GetIdx()) != oxygen_atom.GetIdx() and bond.GetBondTypeAsDouble() == 2.0:
                                                        return True
                                    return False
                                
                                # Check ester linkages
                                if is_ester_linked(sn1_carbon) and is_ester_linked(sn2_carbon):
                                    diacylglycerol_found = True
                                
                                # Check if sn-3 oxygen is connected to phosphate
                                sn3_oxygens = [nbr for nbr in sn3_carbon.GetNeighbors() if nbr.GetAtomicNum() == 8]
                                for oxy_atom in sn3_oxygens:
                                    if oxy_atom.GetIdx() == oxy.GetIdx():
                                        # Oxygen connected to phosphate
                                        # Now check if phosphate is connected to another glycerol (headgroup)
                                        for oxy2 in phosphate_oxygens:
                                            if oxy2.GetIdx() == oxy_atom.GetIdx():
                                                continue
                                            # Check if oxy2 is connected to a glycerol
                                            oxy2_neighbors = [nbr for nbr in oxy2.GetNeighbors() if nbr.GetIdx() != phospho_atom.GetIdx()]
                                            if oxy2_neighbors and oxy2_neighbors[0].GetAtomicNum() == 6:
                                                # Potential headgroup glycerol
                                                # Check for three-carbon chain starting from this carbon
                                                paths2 = Chem.FindAllPathsOfLengthN(mol, 3, useBonds=True, useHs=False)
                                                for path2 in paths2:
                                                    if oxy2_neighbors[0].GetIdx() in path2:
                                                        atoms_in_path2 = [mol.GetAtomWithIdx(idx) for idx in path2]
                                                        if all(atom.GetAtomicNum() == 6 for atom in atoms_in_path2):
                                                            # Check if each carbon has an oxygen attached
                                                            for atom in atoms_in_path2:
                                                                oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                                                                if not oxygens:
                                                                    break  # No oxygen attached
                                                            else:
                                                                headgroup_glycerol_found = True
                                                                break
                                            if headgroup_glycerol_found:
                                                break
                                    if headgroup_glycerol_found:
                                        break
                                if diacylglycerol_found and headgroup_glycerol_found:
                                    return True, "Molecule matches phosphatidylglycerol structure"
        # End of phosphate group iteration

    # If we reach this point, no matching structure was found
    return False, "Molecule does not match phosphatidylglycerol structure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phosphatidylglycerol',
        'definition': 'A glycerophosphoglycerol that is glycerol in which the hydrogen of one of the primary hydroxy groups has been replaced by a phosphatidyl group.',
    },
}