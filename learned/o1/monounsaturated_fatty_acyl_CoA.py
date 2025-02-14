"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety using its SMILES string (simplified for matching)
    coa_smiles = 'CC(C)(COP(=O)(OP(=O)(O)OCC1OC(C(O)C1OP(=O)(O)O)N1C=NC2=C1N=CN=C2N)(O)=O)OC(=O)C(NCCSC)C(=O)NCCCN'
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Invalid CoA SMILES string"
    
    # Check if the molecule contains the CoA moiety
    if not mol.HasSubstructMatch(coa_mol):
        return False, "CoA moiety not found"

    # Define the thioester pattern
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Process each thioester linkage
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[2]

        # Break the bond between carbonyl carbon and sulfur to isolate the acyl chain
        mol_copy = Chem.RWMol(mol)
        bond = mol_copy.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
        if bond is None:
            continue  # Bond not found, skip
        mol_copy.RemoveBond(carbonyl_c_idx, sulfur_idx)
        
        # Get the fragments and the mapping of atom indices
        frags = Chem.GetMolFrags(mol_copy, asMols=False, sanitizeFrags=True)
        # frags is a tuple of tuples containing atom indices for each fragment
        # Identify the fragment containing the carbonyl carbon
        acyl_chain_atoms = None
        for frag in frags:
            if carbonyl_c_idx in frag:
                acyl_chain_atoms = frag
                break

        if acyl_chain_atoms is None:
            continue  # Acyl chain not found, proceed to next thioester linkage

        # Create submol for the acyl chain
        acyl_chain_mol = Chem.RWMol()
        atom_mapping = {}
        for idx in acyl_chain_atoms:
            atom = mol.GetAtomWithIdx(idx)
            new_atom = Chem.Atom(atom.GetAtomicNum())
            acyl_chain_mol_idx = acyl_chain_mol.AddAtom(new_atom)
            atom_mapping[idx] = acyl_chain_mol_idx

        # Add bonds
        for idx in acyl_chain_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx in acyl_chain_atoms and idx < nbr_idx:
                    bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                    acyl_chain_mol.AddBond(atom_mapping[idx], atom_mapping[nbr_idx], bond.GetBondType())
        
        # Count the number of carbon-carbon double bonds in the acyl chain
        double_bonds = 0
        for bond in acyl_chain_mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                # Count only carbon-carbon double bonds
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    double_bonds += 1

        # Check if the acyl chain is predominantly composed of carbons
        c_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
        h_count = sum(1 for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 1)
        heteroatom_count = acyl_chain_mol.GetNumAtoms() - c_count - h_count

        if heteroatom_count > 0:
            return False, "Acyl chain contains heteroatoms, not a fatty acyl chain"

        # Decide based on the number of double bonds
        if double_bonds == 1:
            return True, "Contains one carbon-carbon double bond in the fatty acyl chain"
        elif double_bonds == 0:
            return False, "No carbon-carbon double bonds in the fatty acyl chain"
        else:
            return False, f"Contains {double_bonds} carbon-carbon double bonds in the fatty acyl chain"

    return False, "No suitable thioester linkage associated with CoA found"