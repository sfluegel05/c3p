"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: Glucosinolate
Definition: Water-soluble anionic substituted thioglucosides. 
Glucosinolates have a central C atom which is bonded via an S atom to a glycone group
and via an N atom to a sulfonated oxime group, and which also carries a side-group.
The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.
Note: The anti configuration is not explicitly verified.
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    The molecule is required to have:
      1. A sulfonated oxime group. Our SMARTS pattern (ignoring minor charge variations) is:
             C(=N[O]S(=O)(=O)[O-]?)
         This identifies a C with a double bond to an N, which is connected to an O then to a sulfonyl group.
      2. The central carbon (from above) must have at least three substituents (i.e. including a side chain).
      3. One of the substituents must be a sulfur atom attached to a glycone (sugar) moiety.
         We detect the glycone using a simplified SMARTS for a pyranose ring:
             C1OC(C(O)C(O)C(O)C1O)
         but we perform the matching ignoring stereochemistry.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as glucosinolate, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define the SMARTS for the sulfonated oxime group.
    # This pattern allows for optional charge on the terminal oxygen.
    oxime_smarts = "C(=N[O]S(=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # Unexpected error, should not happen.
    
    # Search for the oxime substructure in the molecule.
    oxime_matches = mol.GetSubstructMatches(oxime_query)
    if not oxime_matches:
        return False, "No sulfonated oxime substructure found."
    
    # Define the SMARTS for a glycone (pyranose sugar) ring.
    # We use a common pyranose representation.
    glycone_smarts = "C1OC(C(O)C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None
    
    # Find glycone matches in the molecule.
    # Here we ignore chirality in matching to allow for different stereo representations.
    glycone_matches = mol.GetSubstructMatches(glycone_query, useChirality=False)
    if not glycone_matches:
        return False, "No glycone (sugar) ring found in the molecule."
    
    # For each detected sulfonated oxime, check that the central carbon:
    #   1. Has at least three substituents (indicating attachment of a side chain).
    #   2. Has a neighbor sulfur atom that connects to the glycone structure.
    found_candidate = False
    reasons = []
    for match in oxime_matches:
        # In our SMARTS the first atom (index 0) is the central carbon.
        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        
        if central_c.GetDegree() < 3:
            reasons.append("Central carbon in sulfonated oxime has fewer than 3 substituents.")
            continue
        
        # Look among neighbors of the central carbon for a sulfur atom.
        sulfur_neighbor = None
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur
                sulfur_neighbor = nbr
                break
        if sulfur_neighbor is None:
            reasons.append("Central carbon not attached to any sulfur atom.")
            continue
        
        # Check if this sulfur is attached (by any bond) to an atom that is part of a glycone match.
        # We iterate through each neighbor of S (excluding the central C already inspected).
        glycone_connected = False
        for s_nbr in sulfur_neighbor.GetNeighbors():
            if s_nbr.GetIdx() == central_c_idx:
                continue
            # Check if s_nbr lies in any of the glycone substructure matches.
            for gly_match in glycone_matches:
                if s_nbr.GetIdx() in gly_match:
                    glycone_connected = True
                    break
            if glycone_connected:
                break
        
        if not glycone_connected:
            reasons.append("Sulfur neighbor is not connected to a glycone ring.")
            continue
        
        # Found a candidate matching the criteria.
        found_candidate = True
        break
    
    if not found_candidate:
        full_reason = " ; ".join(reasons) if reasons else "No valid glucosinolate substructure found."
        return False, full_reason
    
    return True, "Molecule contains a sulfonated oxime group on a central carbon (with a side chain) and a sulfur attached to a glycone (sugar) moiety."

# Example usage:
if __name__ == "__main__":
    # Test a known glucosinolate: 4-(Methylsulfinyl)but-3-enylglucosinolate.
    test_smiles = "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C"
    result, reason = is_glucosinolate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)