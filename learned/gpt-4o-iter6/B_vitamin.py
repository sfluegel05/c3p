"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define enhanced SMARTS patterns for key structural features of B vitamins
    try:
        # Thiamine (B1) patterns including thiazole and pyrimidine
        thiamine_patterns = [
            Chem.MolFromSmarts("Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"),  # Previously used primary pattern
            Chem.MolFromSmarts("C[C@H]1CSC=NC1")  # Thiazole ring
        ]
        # Riboflavin (B2), capturing isoalloxazine ring
        riboflavin_patterns = [
            Chem.MolFromSmarts("Cn1cnc2c1c(=O)n(c(=O)n2C)C")  # Isoalloxazine core simplified
        ]
        # Niacin (B3), nicotinic acid motif
        niacin_patterns = [
            Chem.MolFromSmarts("n1ccc(C(=O)[O,N])cn1")
        ]
        # Pantothenic acid (B5), core structure
        pantothenic_patterns = [
            Chem.MolFromSmarts("OC(=O)C(C(C(=O)O)NC(=O)C(C)O)CO")
        ]
        # Pyridoxine (B6), pyridine derivative
        pyridoxine_patterns = [
            Chem.MolFromSmarts("C1=NC=C(C(=C1CO)CO)CO")
        ]
        # Biotin (B7), characteristic urea and thiophene
        biotin_patterns = [
            Chem.MolFromSmarts("C1(C(=O)NC2C(S1)CCC2)CCC(=O)[O,N]")
        ]
        # Folate (B9), capturing pterin and glutamic acid chains
        folate_patterns = [
            Chem.MolFromSmarts("Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1"),
            Chem.MolFromSmarts("C1=CC(=CC=C1)C(=O)N[C@@H](CCC(=O)[O,N])C(=O)[O,N]")  # Tail structure
        ]
        # Cobalamin (B12), cobalamin core feature, which is extremely complex 
        # Here, we use a basic identifier, real-world application may need specialized checks
        cobalamin_patterns = [
            Chem.MolFromSmarts("[Co]")  # Contains cobalt ion typical of B12
        ]

        vitamin_patterns = {
            "Thiamine (B1)": thiamine_patterns,
            "Riboflavin (B2)": riboflavin_patterns,
            "Niacin (B3)": niacin_patterns,
            "Pantothenic acid (B5)": pantothenic_patterns,
            "Pyridoxine (B6)": pyridoxine_patterns,
            "Biotin (B7)": biotin_patterns,
            "Folate (B9)": folate_patterns,
            "Cobalamin (B12)": cobalamin_patterns
        }

        for vitamin_name, patterns in vitamin_patterns.items():
            if any(mol.HasSubstructMatch(pattern) for pattern in patterns):
                return True, f"Matches {vitamin_name} structure"

        return False, "No match for B vitamin structures found in SMILES"
    except Exception as e:
        return None, f"Error during pattern matching: {str(e)}"