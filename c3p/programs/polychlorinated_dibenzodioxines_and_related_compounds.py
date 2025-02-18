"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds (persistent organic pollutants)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and related compounds.
    These are organochlorine/bromine compounds with dioxin, dibenzofuran, or biphenyl cores,
    where ALL substituents on the core rings are halogens (Cl/Br).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Must contain at least two chlorine or bromine atoms
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    if halogen_count < 2:
        return False, f"Insufficient halogens ({halogen_count} found, minimum 2)"

    # Core structure patterns
    dioxin_core = Chem.MolFromSmarts('O1c2ccccc2Oc2ccccc12')  # Dibenzodioxin
    furan_core = Chem.MolFromSmarts('c1ccc2c(c1)oc1ccccc12')  # Dibenzofuran
    biphenyl_core = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')   # Biphenyl

    # Check each core pattern for valid substitution
    for core in [dioxin_core, furan_core, biphenyl_core]:
        matches = mol.GetSubstructMatches(core)
        for match in matches:
            core_atoms = set(match)
            valid = True
            
            # Check all substituents on core's aromatic atoms are Cl/Br
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                if not atom.GetIsAromatic():
                    continue  # Only check aromatic positions
                    
                for neighbor in atom.GetNeighbors():
                    # Skip core-internal bonds
                    if neighbor.GetIdx() in core_atoms:
                        continue
                        
                    # Substituents must be Cl/Br
                    if neighbor.GetAtomicNum() not in [17, 35]:
                        valid = False
                        break
                
                if not valid:
                    break
            
            if valid:
                core_type = "dioxin" if core == dioxin_core else \
                            "dibenzofuran" if core == furan_core else "biphenyl"
                return True, f"Halogenated {core_type} core with all substituents as Cl/Br"

    return False, "No valid dioxin/dibenzofuran/biphenyl core with halogen-only substituents"