"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones are defined by the presence of a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core,
    which is represented here by a chromanone substructure with an additional aryl substituent.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a flavanone, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern to capture a simplified chromanone (benzopyran-4-one) core.
    # We include an atom mapping on the second carbon (using [C:2]) which is the typical site of aryl substitution.
    core_smarts = "C1[C:2]C(=O)Oc2ccccc12"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern for flavanone core"

    # Search for the core pattern in the molecule.
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Flavanone core structure (chromanone) not found"

    # For at least one core match, check if the atom corresponding to the mapped [C:2] is connected to an additional aromatic substituent.
    # Since the SMARTS uses a mapping we expect the second atom in the query to be at index 1 in each match tuple.
    # (Note: RDKit returns only indices corresponding to the pattern order without the explicit mapping numbers;
    # so we assume that atom at match[1] denotes the flavanone's C-2.)
    for match in matches:
        # Get the index of the candidate C-2 atom in the molecule
        c2_idx = match[1]
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        # Look at neighbors (outside of the core match) to see if an aromatic ring (aryl substituent) is attached
        has_extra_aryl = False
        for neighbor in c2_atom.GetNeighbors():
            # if neighbor is not in the core match, then it could be the extra aryl moiety.
            if neighbor.GetIdx() not in match and neighbor.GetIsAromatic():
                has_extra_aryl = True
                break
        if has_extra_aryl:
            return True, "Molecule contains the flavanone core with an extra aryl substituent (at C-2)"
    return False, "Flavanone core found but missing a 2-aryl substituent"

# Example usage (you can remove or comment this out if using it as a module)
if __name__ == "__main__":
    examples = {
        "eriodictyol": "Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",
        "sakuranetin": "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",
        "hesperetin": "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1",
        "(R)-naringenin": "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    }
    for name, smi in examples.items():
        result, reason = is_flavanones(smi)
        print(f"{name}: {result} ({reason})")