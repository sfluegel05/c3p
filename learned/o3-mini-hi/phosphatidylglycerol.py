"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.
Improved criteria:
  1. Contains exactly one phosphorus atom.
  2. Contains a phosphoglycerol headgroup matching the SMARTS pattern "P(=O)(O)(OCC(O)CO)"
     (this ensures one oxygen substituent extends to a glycerol fragment).
  3. Contains exactly two acyl ester groups (matching "OC(=O)") that are not part of the headgroup.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    A phosphatidylglycerol is defined by having exactly one phosphorous,
    a phosphoglycerol headgroup (in which one primary hydroxyl in glycerol is replaced by a phosphatidyl group),
    and two acyl chains attached as ester groups.
    
    We require:
     - Exactly one phosphorus atom.
     - A unique phosphoglycerol headgroup substructure (SMARTS: "P(=O)(O)(OCC(O)CO)").
     - Exactly two ester groups (SMARTS: "OC(=O)") that are not part of the headgroup.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one phosphorus atom (atomic number 15)
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(phosphorus_atoms)}"
    
    # Define the SMARTS for the phosphoglycerol headgroup.
    # This pattern looks for a phosphorus double-bonded to oxygen, bearing two other oxygens,
    # one of which is connected to a 2-carbon chain with two hydroxy groups: OCC(O)CO.
    pg_head_smarts = "P(=O)(O)(OCC(O)CO)"
    pg_head = Chem.MolFromSmarts(pg_head_smarts)
    pg_matches = mol.GetSubstructMatches(pg_head, useChirality=False)
    if len(pg_matches) != 1:
        return False, f"Phosphoglycerol headgroup pattern not found exactly once (found {len(pg_matches)} occurrence(s))"
    
    # Gather indices of atoms that are part of the headgroup match
    pg_atom_indices = set(pg_matches[0])
    
    # Define the SMARTS for acyl ester groups
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Filter ester matches to those that are NOT part of the headgroup match.
    valid_ester_matches = []
    for match in ester_matches:
        # If none of the atoms in the ester match is in the headgroup, then count this as an acyl ester.
        if not set(match) & pg_atom_indices:
            valid_ester_matches.append(match)
    
    if len(valid_ester_matches) != 2:
        return False, f"Expected 2 acyl ester groups (outside headgroup), found {len(valid_ester_matches)}"
    
    # Optional additional checks could be performed here (e.g. verifying the length of acyl chains)
    # For now, if our headgroup and two acyl esters are found, we classify as phosphatidylglycerol.
    
    return True, "Molecule contains a phosphoglycerol headgroup and exactly two acyl ester chains consistent with phosphatidylglycerol"

# Example usage:
if __name__ == "__main__":
    # Test with one example: PG(8:0/8:0)
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)