"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57844 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is an amino acid zwitterion obtained by transfer of a proton from 
    the carboxy group to the amino group of any alpha-amino acid; major species at physiological pH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms to find protonated amino groups [NH3+]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Check if nitrogen has 4 bonds (including implicit hydrogens)
            if atom.GetTotalDegree() == 4:
                # Should have 3 hydrogens
                if atom.GetTotalNumHs(includeNeighbors=True) == 3:
                    # Get alpha carbon connected to this nitrogen
                    neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
                    for alpha_carbon in neighbors:
                        # Check if alpha carbon is sp3 hybridized
                        if alpha_carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                            continue
                        # Check if alpha carbon is connected to carboxylate group [C(=O)[O-]]
                        has_carboxylate = False
                        for nbr2 in alpha_carbon.GetNeighbors():
                            if nbr2.GetAtomicNum() == 6 and nbr2.GetIdx() != atom.GetIdx():
                                # Check for carboxylate carbon
                                carboxylate_carbon = nbr2
                                oxygens = [o for o in carboxylate_carbon.GetNeighbors() if o.GetAtomicNum() == 8]
                                if len(oxygens) == 2:
                                    charges = [o.GetFormalCharge() for o in oxygens]
                                    # One oxygen must have a negative charge
                                    if -1 in charges:
                                        # Carboxylate group found
                                        has_carboxylate = True
                                        break
                        if has_carboxylate:
                            # Check for side chain (any other group connected to alpha carbon)
                            side_chain_atoms = [nbr for nbr in alpha_carbon.GetNeighbors()
                                                if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() != 6]
                            if len(side_chain_atoms) >= 1:
                                return True, "Contains alpha-amino-acid zwitterion structure"
                            else:
                                # For glycine, side chain is hydrogen
                                h_count = sum(1 for nbr in alpha_carbon.GetNeighbors() if nbr.GetAtomicNum() == 1)
                                if h_count >= 1:
                                    return True, "Contains alpha-amino-acid zwitterion structure (glycine)"
    return False, "Alpha-amino-acid zwitterion pattern not found"