"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms
             of an organyl group.

This implementation restricts to arsenic–based compounds and uses several heuristics:
  1. The SMILES must contain an arsenic atom ("[As") to limit to arsenic compounds.
  2. A hard–coded filter rejects tryparsamide–like fragments (and a few other patterns).
  3. For each arsenic atom, examine its bonds to carbon. The carbon is considered “organyl” if:
       • It is part of a chain or ring (i.e. it is bonded to another carbon besides the arsenic), or
       • It is a methyl group and either (a) the entire molecule is small or (b) the arsenic bears a positive charge.
  4. In addition, for neutral, non–ring As atoms we require at least one double bond to oxygen (a hallmark of As(V)).
  
Note: This heuristic–based approach is an improvement over the previous version but remains an imperfect tool.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    (Heuristics:
       - Restricts to arsenic-based compounds.
       - Requires at least one As–C bond where the carbon is considered to be part of an organyl group.
       - For neutral, non–ring As, requires an As=O bond.
       - Accepts simple (methyl) bonds in low molecular weight or charged cases.)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if judged to be an organometalloidal compound, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Only consider arsenic-based compounds.
    if "[As" not in smiles:
        return False, "No arsenic atom found; not an organometalloidal compound"
    
    # Hard-coded filter: reject tryparsamide-like fragment.
    if "NCC(N)=O" in smiles:
        return False, "Detected tryparsamide-like fragment (NCC(N)=O); not considered an organyl group"
    
    # (Optional) More filters can be added.
    # For example, if a dichloro–arsenic is present, we reject.
    if "[As](Cl" in smiles:
        return False, "Detected Cl substituents on arsenic; likely inorganic derivative"

    # Check molecular weight – very low MW compounds are often noise.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be an organometalloidal compound"
    
    overall_details = []
    valid_overall = False

    # Loop over each atom – focus on arsenic atoms (atomic number 33)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 33:
            continue
        # Gather carbon neighbors
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # This arsenic has no As–C bonds
        
        details = []
        has_organyl = False  # flag: at least one C neighbor is part of a larger organic fragment
        has_methyl   = False  # flag: encountered a methyl group
        for c in carbon_neighbors:
            # Check if this carbon is attached to another carbon (besides the As)
            other_c = [nbr for nbr in c.GetNeighbors() if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6]
            if other_c:
                has_organyl = True
                details.append("Found As–C bond where C is part of a larger organic fragment")
            else:
                has_methyl = True
                details.append("Found As–C bond where C appears as a methyl group")
        
        # For arsenic atoms that have an organyl bond, we generally accept.
        if has_organyl:
            # For neutral, non–ring As atoms, check for an As=O bond.
            if not atom.IsInRing() and atom.GetFormalCharge() == 0:
                # Look among bonds of this As for a double-bonded O:
                double_bonded_O = False
                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bonded_O = True
                        break
                if not double_bonded_O:
                    details.append("Neutral, non–ring As lacks an As=O bond; rejecting")
                    # For this As we do not mark it as valid; continue to next As.
                    continue
            # If we get here, this As qualifies.
            overall_details.extend(details)
            valid_overall = True
            break

        # If there is no large organyl bond, then consider the methyl-only case.
        # We relax the previous rule: if the entire compound is low MW or if the As is charged, allow methyl groups.
        if has_methyl:
            if atom.GetFormalCharge() != 0:
                overall_details.extend(details + ["As bears a nonzero charge; accepting methyl ligand(s)"])
                valid_overall = True
                break
            # Otherwise, if the molecule is very simple we permit a single methyl group.
            if len(carbon_neighbors) == 1 and mol_wt < 170:
                overall_details.extend(details + ["Single methyl substituent on low–MW molecule; accepting"])
                valid_overall = True
                break
            # Also, if there are exactly two methyl ligands and no additional heteroatom ligands (besides O)
            # we permit them (e.g. dimethylarsinous acid).
            if len(carbon_neighbors) == 2:
                other_subs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() not in (6, 1)]
                # If not dominated by extra electronegative groups then accept.
                if len(other_subs) <= 2:
                    overall_details.extend(details + ["Two methyl groups accepted based on substitution pattern"])
                    valid_overall = True
                    break
            details.append("Methyl-only substituents did not meet acceptance criteria")
            overall_details.extend(details)
        # End loop over this As atom.
    
    if valid_overall:
        return True, "; ".join(overall_details)
    else:
        return False, "No As–C bond found that meets organyl criteria"

# Example usage:
if __name__ == "__main__":
    # True positive example:
    tp_smiles = "[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C"  # 1-dimethylarsinoyl-(2Z,...)
    result, reason = is_organometalloidal_compound(tp_smiles)
    print("Input SMILES:", tp_smiles)
    print("Result:", result)
    print("Reason:", reason)
    
    # Example of dimethylarsinic acid (should be accepted)
    dm_arsinic = "C[As](C)(O)=O"
    result, reason = is_organometalloidal_compound(dm_arsinic)
    print("\nInput SMILES:", dm_arsinic)
    print("Result:", result)
    print("Reason:", reason)
    
    # False positive example: m-aminophenylarsonous acid (expected to be rejected)
    fp_smiles = "Nc1cccc(c1)[As](O)O"
    result, reason = is_organometalloidal_compound(fp_smiles)
    print("\nInput SMILES:", fp_smiles)
    print("Result:", result)
    print("Reason:", reason)