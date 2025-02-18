"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'
This function checks if the molecule (given as a SMILES string) contains at least one –OH 
group on an sp3, non-aromatic carbon whose immediate carbon neighbors (aside from the –OH)
are not aromatic. This extra check helps to avoid classifying multi-functional or aromatic
compounds that merely contain an aliphatic -OH group.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined here as one containing at least one hydroxyl (-OH) group 
    attached to a saturated (sp3), non-aromatic carbon, whose local environment (other immediate
    heavy-atom links to that carbon) is not directly connected to aromatic systems.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Iterate over oxygen atoms that might be part of hydroxyl groups.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Check that the oxygen has at least one hydrogen attached 
            # (could be explicit or implicit)
            if atom.GetTotalNumHs() < 1:
                continue  # not an -OH group, skip
            
            # For every neighbor of the oxygen, if it is a carbon then it is the candidate atom.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    # First, ensure that the carbon is sp3 and non-aromatic.
                    if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or nbr.GetIsAromatic():
                        continue  # not an aliphatic carbon
                    
                    # Now check the immediate environment of this candidate carbon.
                    # We allow the oxygen (our -OH) but check that no other attached heavy atom 
                    # (non-hydrogen) is aromatic.
                    passes_environment = True
                    for sub_nbr in nbr.GetNeighbors():
                        if sub_nbr.GetIdx() == atom.GetIdx():
                            continue  # skip our -OH oxygen
                        # Only check non-hydrogen atoms.
                        if sub_nbr.GetAtomicNum() != 1:
                            # If any other neighbor is part of an aromatic ring, discard candidate.
                            if sub_nbr.GetIsAromatic():
                                passes_environment = False
                                break
                    if passes_environment:
                        return True, "Found -OH group attached to an aliphatic (sp3, non-aromatic) carbon with a non-aromatic local environment"
    
    return False, "No qualifying aliphatic -OH group found; either there is no hydroxyl group or its environment is too decorated by aromatic/heteroatom functionality"


# (Optional) Testing examples – these are not required in the final program.
if __name__ == "__main__":
    test_smiles = [
        "O=C1OC([C@@H](O)\\C=C/C=C/C)CC1",  # Sapinofuranone A (should be True)
        "CCCCCCC(C)O",                     # octan-2-ol (should be True)
        "O1C(C(OCC=C(C)C)COC2=C1C=C(C=C2)\\C=C\\C=C(/CO)C(\\C(OC)=O)=C\\OC)(C)C",  # Hydroxystrobilurin D (True)
        "OC1=CC=CC=C1",                    # Phenol: -OH on aromatic C (should be False)
    ]
    
    for sm in test_smiles:
        is_valid, reason = is_aliphatic_alcohol(sm)
        print(f"SMILES: {sm}\nClassification: {is_valid}\nReason: {reason}\n")