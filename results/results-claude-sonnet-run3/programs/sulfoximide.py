from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfoximide(smiles: str):
    """
    Determines if a molecule contains a sulfoximide group (RS(=O)=NR).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains sulfoximide group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for sulfoximide group: RS(=O)=NR
    # [S;X4] - sulfur with 4 bonds
    # (=O) - double bonded oxygen
    # (=N) - double bonded nitrogen or single bonded nitrogen (since SMILES can represent it either way)
    sulfoximide_pattern = Chem.MolFromSmarts('[#6,#1][S;X4](=[O,N])([#6,#1])=[O,N]')
    
    if mol.HasSubstructMatch(sulfoximide_pattern):
        matches = mol.GetSubstructMatches(sulfoximide_pattern)
        
        # Verify that both O and N are present in the matched structure
        for match in matches:
            s_atom = mol.GetAtomWithIdx(match[1])  # Sulfur atom
            s_neighbors = s_atom.GetNeighbors()
            
            has_oxygen = False
            has_nitrogen = False
            
            for neighbor in s_neighbors:
                if neighbor.GetSymbol() == 'O':
                    has_oxygen = True
                elif neighbor.GetSymbol() == 'N':
                    has_nitrogen = True
                    
            if has_oxygen and has_nitrogen:
                # Get R groups
                r_groups = []
                for neighbor in s_neighbors:
                    if neighbor.GetSymbol() not in ['O', 'N']:
                        r_groups.append(neighbor.GetSymbol())
                        
                return True, f"Contains sulfoximide group with R groups: {', '.join(set(r_groups))}"
    
    return False, "No sulfoximide group found"
# Pr=1.0
# Recall=1.0