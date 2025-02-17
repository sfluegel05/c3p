"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: Glucosinolate 
Definition: Water-soluble anionic substituted thioglucosides. 
Glucosinolates have a central C atom which is bonded via an S atom to a glycone group 
and via an N atom to a sulfonated oxime group, and which also carries a side-group.
The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.
Note: The stereochemical “anti” configuration is not explicitly verified.
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    We require that the molecule contains:
    1. A sulfonated oxime group attached to a central carbon. This is approximated with the SMARTS:
         C(=N[O][S](=O)(=O)[O-]?)
       which looks for a C double-bonded to an N which is bound to an O then to a sulfonyl group.
    2. The same central carbon must also be attached to a sulfur atom that is part of a glycone (sugar) ring.
       We use a simplified sugar SMARTS for a pyranose ring (typical for glucose):
         C1OC(C(O)C(O)C(O)C1O)
    
    It is also verified that the central carbon has at least three substituents (indicating that it 
    also carries a side-chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as glucosinolate, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define the SMARTS for the sulfonated oxime group.
    # This pattern looks for a carbon atom double-bonded to an N,
    # with that N attached to an oxygen then to a sulfonyl group (OS(=O)(=O)[O-] or neutral OS(O)(=O)=O)
    oxime_smarts = "C(=N[O][S](=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # Should not happen if SMARTS is valid.
    
    # Search for the sulfonated oxime substructure.
    oxime_matches = mol.GetSubstructMatches(oxime_query)
    if not oxime_matches:
        return False, "No sulfonated oxime substructure found."
    
    # Define the SMARTS for a glycone (sugar) ring.
    # This simplified pyranose pattern covers many typical glucosides.
    glycone_smarts = "C1OC(C(O)C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None
    
    # For each sulfonated oxime match, check that:
    # 1. The central carbon (atom 0 in our oxime SMARTS) has at least three substituents.
    # 2. One of the neighbors of this carbon is a sulfur atom that is part of a glycone substructure.
    found_candidate = False
    reason_details = []
    for match in oxime_matches:
        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        # Check that central C has at least 3 neighbors (S, N and a side-group)
        if central_c.GetDegree() < 3:
            reason_details.append("Central C in oxime does not have at least 3 substituents.")
            continue  # Try next match
        
        # Look through neighbors of the central carbon for an S atom.
        neighbor_s_found = False
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # 16 is sulfur
                # Check if the candidate S is part of a glycone ring.
                # We get all matches of the glycone pattern in the molecule
                glycone_matches = mol.GetSubstructMatches(glycone_query)
                for gly_match in glycone_matches:
                    if nbr.GetIdx() in gly_match:
                        neighbor_s_found = True
                        break
            if neighbor_s_found:
                break
        
        if not neighbor_s_found:
            reason_details.append("Central C not attached to a sulfur that is part of a glycone ring in this match.")
            continue
        
        # If we reach here, we have found an oxime match with the correct central carbon that has:
        # -- at least 3 substituents,
        # -- an S neighbor that appears in a sugar ring.
        found_candidate = True
        break
    
    if not found_candidate:
        # If we did not find any candidate that satisfies both requirements, compile a reason.
        full_reason = " ; ".join(reason_details) if reason_details else "No valid glucosinolate substructure found."
        return False, full_reason
    
    # (Optional further checks could include molecular weight or count of OH groups for water solubility.)
    # Note: We are not verifying the anti stereochemical configuration across the C=N double bond.
    
    return True, "Molecule contains a sulfonated oxime group attached to a central C that is also linked to a thioglucoside (glycone) group with a side-chain."

# Example usage:
if __name__ == "__main__":
    # Test one known glucosinolate SMILES: 4-(Methylsulfinyl)but-3-enylglucosinolate
    test_smiles = "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C"
    result, reason = is_glucosinolate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)