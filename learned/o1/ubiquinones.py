"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone is defined as a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core SMILES for 2,3-dimethoxy-5-methylbenzoquinone
    core_smiles = 'COC1=C(C=C(C(=O)C1=O)C)OC'
    core_mol = Chem.MolFromSmiles(core_smiles)
    if core_mol is None:
        return False, "Unable to create core molecule"

    # Check for the core pattern in the molecule
    matches = mol.GetSubstructMatches(core_mol)
    if not matches:
        return False, "Core 2,3-dimethoxy-5-methylbenzoquinone not found"

    # For each match, find if there is a side chain attached to the core at position 6
    for match in matches:
        matched_atom_ids = set(match)
        # Iterate over the atoms in the core and look for attachments
        for core_atom_idx, mol_atom_idx in enumerate(match):
            mol_atom = mol.GetAtomWithIdx(mol_atom_idx)
            # Check for atoms connected to the core atom that are not in the core
            for bond in mol_atom.GetBonds():
                neighbor = bond.GetOtherAtom(mol_atom)
                if neighbor.GetIdx() not in matched_atom_ids:
                    # Found an external attachment, could be the side chain
                    side_chain_atoms = set()
                    atoms_to_check = [neighbor]
                    while atoms_to_check:
                        current_atom = atoms_to_check.pop()
                        if current_atom.GetIdx() in side_chain_atoms:
                            continue
                        side_chain_atoms.add(current_atom.GetIdx())
                        for bond2 in current_atom.GetBonds():
                            next_atom = bond2.GetOtherAtom(current_atom)
                            if next_atom.GetIdx() not in side_chain_atoms and next_atom.GetIdx() not in matched_atom_ids:
                                atoms_to_check.append(next_atom)
                    # Analyze the side chain
                    if len(side_chain_atoms) == 0:
                        continue  # Empty side chain
                    side_chain_mol = Chem.PathToSubmol(mol, list(side_chain_atoms))
                    if side_chain_mol is None:
                        continue  # Unable to isolate side chain
                    # Check if the side chain is a polyprenoid
                    num_carbons = sum(1 for atom in side_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
                    if num_carbons < 5:
                        continue  # Side chain too short
                    num_double_bonds = sum(1 for bond in side_chain_mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
                    if num_double_bonds < 1:
                        continue  # Side chain lacks double bonds
                    # Check for isoprene units
                    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
                    isoprene_matches = side_chain_mol.GetSubstructMatches(isoprene_unit)
                    if len(isoprene_matches) < 1:
                        continue  # Side chain does not contain isoprene units
                    # If all checks passed, return True
                    return True, "Contains 2,3-dimethoxy-5-methylbenzoquinone core with polyprenoid side chain at position 6"

    # If no matches with side chain found, return False
    return False, "Core found but no suitable side chain attached at position 6"

__metadata__ = {   
    'chemical_class': {   
        'name': 'ubiquinones',
        'definition': 'Any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone; one of a group of naturally occurring homologues. The redox-active quinoid moiety usually carries a polyprenoid side chain at position 6, the number of isoprenoid units in which is species-specific. Ubiquinones are involved in the control of mitochondrial electron transport, and are also potent anti-oxidants.',
    }
}