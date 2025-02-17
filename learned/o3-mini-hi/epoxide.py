"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
In a valid epoxide, the ring contains exactly one oxygen and two carbons, 
all in a 3‐membered ring with single bonds and sp³ hybridization.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxide(smiles: str):
    """
    Determines if the molecule represented by the SMILES string is an epoxide.
    An epoxide is defined as a cyclic ether with a 3-membered ring (one oxygen, two carbons),
    where the atoms in that ring are typically sp³-hybridized and the bonds are single bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an epoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a 3-membered cyclic ether.
    # The pattern "[#6;r3]-[#8;r3]-[#6;r3]" ensures that each atom in the match (two carbons and one oxygen)
    # is in a 3-membered ring.
    epoxide_smarts = "[#6;r3]-[#8;r3]-[#6;r3]"
    epoxide_query = Chem.MolFromSmarts(epoxide_smarts)
    if epoxide_query is None:
        return False, "Failed to create an epoxide SMARTS pattern"
    
    # Find all substructure matches that meet the basic criteria.
    matches = mol.GetSubstructMatches(epoxide_query, useChirality=True)
    valid_count = 0

    # Loop through each match and further validate.
    for match in matches:
        # Retrieve the atoms corresponding to this match.
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Check that each atom is sp3 hybridized (typical for an epoxide).
        if not all(atom.GetHybridization() == rdchem.HybridizationType.SP3 for atom in atoms):
            continue  # Skip if any atom is not sp3
        
        # Check that the three atoms form a closed ring having only single bonds.
        # We examine the bonds between atom0-atom1, atom1-atom2, and atom2-atom0.
        bond1 = mol.GetBondBetweenAtoms(match[0], match[1])
        bond2 = mol.GetBondBetweenAtoms(match[1], match[2])
        bond3 = mol.GetBondBetweenAtoms(match[2], match[0])
        if bond1 is None or bond2 is None or bond3 is None:
            continue
        if (bond1.GetBondType() != rdchem.BondType.SINGLE or 
            bond2.GetBondType() != rdchem.BondType.SINGLE or 
            bond3.GetBondType() != rdchem.BondType.SINGLE):
            continue
        
        # If all criteria are met, count this as a valid epoxide ring.
        valid_count += 1

    if valid_count > 0:
        return True, f"Found epoxide ring(s): {valid_count} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No valid epoxide ring (3-membered cyclic ether with one oxygen and two carbons) found"