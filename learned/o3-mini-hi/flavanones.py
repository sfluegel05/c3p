"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
This pattern is detected by searching for a 3,4-dihydrochroman-4-one core with an extra 
aromatic (aryl) substituent attached at the C-2 position.
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    
    Flavanones are defined by the presence of the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.
    In our approach we search for a heterocycle that:
      - Contains one oxygen (position 1 of the pyran ring)
      - Has two consecutive saturated carbons (positions 2 and 3),
        with the position 2 carbon bearing an extra (non-core) aromatic substituent.
      - Has a ketone (C(=O)) at position 4.
      - Is fused to an aromatic ring (positions 5–9 in the bicyclic core).
    
    The SMARTS below attempts to capture that core:
        O1[C;H]([*])[C;H]C(=O)c2ccccc12
    where the wildcard attached to the C at position 2 is later probed to confirm it is
    connected to an extra aromatic substituent.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 3,4-dihydrochroman-4-one core.
    # Breakdown of the SMARTS:
    #  "O1"       : oxygen atom that starts the ring (position 1)
    #  "[C;H]([*])" : saturated carbon (position 2) bearing some substituent (*)
    #  "[C;H]"    : saturated carbon (position 3)
    #  "C(=O)"    : carbonyl carbon (position 4)
    #  "c2ccccc2" : aromatic ring (the fused A ring)
    #  "12"       : ring closures connecting atom 1 and aromatic ring.
    core_smarts = "O1[C;H]([*])[C;H]C(=O)c2ccccc12"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern for flavanone core"
    
    # Search for the flavanone core in the molecule.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Flavanone core structure (3,4-dihydro-2-aryl-1-benzopyran-4-one) not found"
    
    # For at least one match, check that the C-2 atom (the second atom in the SMARTS match)
    # is attached to an extra aryl (aromatic) substituent not part of the matched core.
    # RDKit returns a tuple of atom indices in the order that the SMARTS pattern defines, so
    # we assume that index 1 corresponds to the C-2 atom.
    for match in core_matches:
        c2_idx = match[1]
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        found_aryl = False
        # Check all neighbors of C-2.
        for neighbor in c2_atom.GetNeighbors():
            # If a neighbor is not part of the core match and is aromatic, it is likely the extra aryl.
            if neighbor.GetIdx() not in match and neighbor.GetIsAromatic():
                found_aryl = True
                break
        if found_aryl:
            return True, "Molecule contains the flavanone core with an extra 2-aryl substituent"
    
    return False, "Flavanone core found but missing a 2-aryl substituent"

# Example usage (can be removed if used as a module)
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