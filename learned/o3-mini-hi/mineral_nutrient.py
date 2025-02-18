"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral nutrient is defined as a mineral (from a set of allowed metals)
that is an inorganic nutrient essential for metabolic/structural function.
Our heuristic idea is three‐part:
  (1) There must be at least one acceptable metal ion. We accept only those
      metal‐containing fragments that are clearly inorganic—in particular,
      either an isolated metal ion (a one‐atom fragment with a positive formal charge)
      or a fragment containing no carbon atoms.
  (2) Any fragment that does not contain an acceptable metal (the counterion)
      must be checked. If the fragment is organic (contains carbon) then it is allowed only if:
         – it is very simple (i.e. a carboxylate with 1–3 carbons) or (if larger) very lipophilic 
           (≥8 carbons) and shows no extra organic (aromatic or non‐carbonyl C=C) features.
      Additionally, we require a slight metal‐dependency: for alkali metals (Na, K, Cs)
      only the simplest (one‐carbon, i.e. formate) organic counterion is tolerated.
  (3) If any non‐water fragment looks overall “organic” (an aromatic or non‐carbonyl double bond between C)
      outside the allowed counterion patterns then we disqualify the entire molecule.
Because the boundaries are not sharp, this is just one possible heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    
    Returns:
        (bool, str): Tuple of classification (True for a mineral nutrient) and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # For assessing small fragments we add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Allowed metal symbols (nutrient metals). Note that later we distinguish alkali metals.
    allowed_metals = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    alkali_metals = {"Na", "K", "Cs"}
    
    # ----------------------------- Helper functions -----------------------------
    def is_water(frag):
        # Identify water: exactly one O and two H (order insensitive)
        atoms = [atom.GetSymbol() for atom in frag.GetAtoms()]
        return len(atoms)==3 and sorted(atoms)==["H", "H", "O"]
    
    def carbon_count(frag):
        return sum(1 for atom in frag.GetAtoms() if atom.GetSymbol() == "C")
    
    def fragment_is_organic(frag):
        """
        Returns True if the fragment has carbon and either contains any aromatic atom,
        or has a double bond between carbons that is not part of a carbonyl.
        """
        if not any(atom.GetSymbol() == "C" for atom in frag.GetAtoms()):
            return False
        # Check for aromatic atoms:
        if any(atom.GetIsAromatic() for atom in frag.GetAtoms()):
            return True
        # Check bonds: if a double bond exists between carbons and is not a carbonyl bond.
        for bond in frag.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetSymbol()=="C" and a2.GetSymbol()=="C":
                    # verify if either is in a carbonyl
                    def in_carbonyl(atom):
                        for nbr in atom.GetNeighbors():
                            b = frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if nbr.GetSymbol() == "O" and b and b.GetBondType() == Chem.BondType.DOUBLE:
                                return True
                        return False
                    if not (in_carbonyl(a1) or in_carbonyl(a2)):
                        return True
        return False

    def is_carboxylate(frag):
        """
        Checks if the fragment contains a carboxylate substructure.
        For our purposes we use the SMARTS: [CX3](=O)[O-] 
        (that is: a trigonal carbon bound to a doubly bonded O and an O-).
        """
        carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O-]")
        return frag.HasSubstructMatch(carboxylate_smarts)
    
    def allowed_organic_counterion(frag, metal_is_alkali: bool):
        """
        Decides if an organic (carbon‐containing) counterion fragment is allowed.
         – It must be very simple (i.e. a carboxylate with only a few carbons)
           OR if it is large (>=8 carbons) then it must not show additional organic features.
         – For alkali metals, only one–carbon carboxylates (i.e. formate) are allowed.
        """
        cnt = carbon_count(frag)
        # First, require that it looks like a carboxylate.
        if not is_carboxylate(frag):
            # If not carboxylate then if it is lipophilic and shows no extra organic features, allow it.
            if cnt >= 8 and not fragment_is_organic(frag):
                return True
            return False
        # It is a carboxylate.
        if metal_is_alkali:
            # For alkali metals, only allow formate (1 carbon).
            return cnt == 1
        else:
            # For non-alkali allowed metals, allow 1-3 carbons.
            return cnt in {1, 2, 3}
    
    # ------------------------- Main classification -------------------------
    
    # Split into fragments (on '.')
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    accepted_metal_fragments = []  # list of (metal_symbol, description) for accepted metal ions
    # We'll also record metal fragments that are found but that may be “borderline”
    metal_reason_msgs = []
    
    # (A) Process each fragment to identify acceptable metal-containing fragments.
    for frag in frags:
        if is_water(frag):
            continue
        frag_atoms = frag.GetAtoms()
        frag_has_allowed_metal = False
        # For each atom in this fragment, if it is an allowed metal…
        for atom in frag_atoms:
            sym = atom.GetSymbol()
            if sym not in allowed_metals:
                continue
            # If the fragment has any carbon, we do not accept the metal in that fragment (it is too “organic”)
            if carbon_count(frag) > 0:
                continue
            # Case 1: The fragment is only one atom and its formal charge is positive.
            if frag.GetNumAtoms() == 1 and atom.GetFormalCharge() > 0:
                accepted_metal_fragments.append((sym, f"Found isolated metal ion: {sym}{atom.GetFormalCharge()}"))
                frag_has_allowed_metal = True
                break
            # Case 2: Multi-atom inorganic fragment (no carbon). For our heuristic we require that at least one metal
            # atom in the fragment has a nonzero formal charge (or the overall fragment charge is positive).
            overall_frag_charge = Chem.GetFormalCharge(frag)
            if overall_frag_charge > 0:
                accepted_metal_fragments.append((sym, f"Found inorganic metal in non-organic fragment: {sym}{'' if atom.GetFormalCharge()==0 else atom.GetFormalCharge()}"))
                frag_has_allowed_metal = True
                break
        if frag_has_allowed_metal:
            # also record a message for this fragment
            pass
    # If no accepted metal-containing fragment is found, disqualify.
    if not accepted_metal_fragments:
        return False, "No acceptable metal nutrient ion found"
    
    # Determine if any accepted metal is alkali (we require a conservative counterion check for them).
    metal_is_alkali = any(m[0] in alkali_metals for m in accepted_metal_fragments)
    
    # (B) Now check counterion (non–metal) fragments.
    counterion_msgs = []
    for frag in frags:
        if is_water(frag):
            continue
        # Skip fragments that contain any allowed metal (they were already processed)
        if any(atom.GetSymbol() in allowed_metals and carbon_count(frag)==0 for atom in frag.GetAtoms()):
            continue
        # For fragments that contain no allowed metal:
        # If the fragment has carbon we check if it is acceptable.
        if carbon_count(frag) > 0:
            # if the fragment appears overtly organic (aromatic or non-carbonyl C=C), disqualify.
            if fragment_is_organic(frag):
                return False, "Disqualifying organic fragment detected: " + Chem.MolToSmiles(frag)
            # Otherwise, if it’s organic but not overtly so, require that it be a simple carboxylate OR big and lipophilic.
            if not allowed_organic_counterion(frag, metal_is_alkali):
                counterion_msgs.append("Disallowed counterion fragment: " + Chem.MolToSmiles(frag))
        # Fragments with no carbon (e.g. chloride, phosphate, etc.) are allowed.
    
    # (C) Summarize metal findings.
    for m in accepted_metal_fragments:
        metal_reason_msgs.append(m[1])
    
    # (D) Final decision: if any counterion issues were found for a molecule with free metal ions then disqualify.
    if counterion_msgs:
        return False, ("Found metal nutrient ion(s): " + " ; ".join(metal_reason_msgs)
                       + " but with counterion issue(s): " + " ; ".join(counterion_msgs))
    
    return True, "Found metal nutrient ion(s): " + " ; ".join(metal_reason_msgs)


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = {
      "Potassium hexachloropalladate(IV)": "[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]",
      "iron(3+) phosphate": "[Fe+3].[O-]P([O-])(=O)[O-]",
      "caesium formate": "[Cs+].[H]C([O-])=O",
      "calcium silicate": "[Ca++].[Ca++].[O-][Si]([O-])([O-])[O-]",
      "zinc nitrate": "[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
      "aluminium sulfate octadecahydrate": "O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
      "magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
      "antimony pentafluoride": "[Sb](F)(F)(F)(F)F",
      "caesium chloride": "[Cl-].[Cs+]",
      "Calcium hypochlorite": "Cl[O-].[Ca+2].Cl[O-]",
      "Lanthanum trichloride": "Cl[La](Cl)Cl",
      "magnesium phosphate": "P([O-])([O-])([O-])=O.P([O-])([O-])([O-])=O.[Mg+2].[Mg+2].[Mg+2]",
      "calcium dichloride": "[Cl-].[Cl-].[Ca++]",
      "trisodium phosphate": "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O",
      "barium sulfate": "[O-]S([O-])(=O)=O.[Ba+2]",
      "aluminium sulfate (anhydrous)": "[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
      "barium nitrate": "[Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
      "calcium monohydroxide": "O[Ca]",
      "calcium sulfate": "[Ca++].[O-]S([O-])(=O)=O",
      "magnesium dipropionate": "[Mg++].CCC([O-])=O.CCC([O-])=O",
      "tricalcium bis(phosphate)": "[Ca++].[Ca++].[Ca++].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
      "disodium hydrogenphosphate": "[Na+].[Na+].OP([O-])([O-])=O",
      "potassium chloride": "[Cl-].[K+]",
      "calcium dihydroxide": "[OH-].[OH-].[Ca++]",
      "magnesium dichloride": "[Mg++].[Cl-].[Cl-]",
      "magnesium sulfate": "[Mg++].[O-]S([O-])(=O)=O",
      "barium acetate": "[Ba++].CC([O-])=O.CC([O-])=O",
      "barium carbonate": "[Ba++].[O-]C([O-])=O",
      "calcium carbonate": "[Ca+2].C(=O)([O-])[O-]",
      "calcium difluoride": "[F-].[F-].[Ca++]",
      "potassium sulfate": "[K+].[K+].[O-]S([O-])(=O)=O",
      "calcium hydrogenphosphate": "[Ca++].[H]OP([O-])([O-])=O",
      # Some false positives (expected to be rejected)
      "antimonic acid": "[H]O[Sb](=O)(O[H])O[H]",
      "aluminium trichloride hexahydrate": "[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].Cl[Al](Cl)Cl",
      "hexaaquamagnesium(2+)": "[H][O]([H])[Mg++]([O]([H])[H])([O]([H])[H])([O]([H])[H])([O]([H])[H])[O]([H])[H]",
      "aluminium hydroxide": "[H]O[Al](O[H])O[H]",
      "sodium hexacyanoferrate(4-)": "[Na+].[Na+].[Na+].[Na+].N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N",
      "Taurocholic acid sodium salt hydrate": "S([O-])(=O)(=O)CCNC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)CC4)[H])C)(C[C@@H]2O)[H])[H])(CC1)[H])C)[H])C.[Na+].O",
      "disodium (alpha-D-galactopyranosyluronate)-(1->4)-alpha-D-galactopyranuronate": "[Na+].[Na+].O[C@H]1O[C@@H]([C@H](O[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)C([O-])=O)[C@H](O)[C@H]1O)C([O-])=O",
      "sodium chlorite": "[Na+].[O-][Cl]=O",
      "sodium hydroxide": "[OH-].[Na+]",
      "sodium 3-aminopropyl 2-acetamido-2-deoxy-alpha-D-glucose-1-phosphate": "[Na+].CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1OP([O-])(=O)OCCCN",
      "potassium tetrabromoaurate": "[K+].Br[Au-](Br)(Br)Br",
      "barium monohydroxide": "O[Ba]",
      "calcium titanate": "[Ca+2].[Ti+4].[O-2].[O-2].[O-2]",
      "aluminium phosphide": "[Al+3].[P-3]",
      "Sodium;(2S)-1-hydroxy-2-[[(2S)-4-methyl-2-[[1-[(2-methylpropan-2-yl)oxycarbonyl]piperidin-4-yl]oxycarbonylamino]pentanoyl]amino]-3-[(3S)-2-oxopyrrolidin-3-yl]propane-1-sulfonate": "[Na+].O=S(=O)(C([C@H](C[C@H]1C(NCC1)=O)NC(C)=O)[O-]",
      "sodium dihydrogenphosphate monohydrate": "P(=O)(O)([O-])O.[Na+].O",
      "sodium glycocholate": "[Na+].[Na+].[O-][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCC(=O)NCC([O-])=O",
      "sodium selenate": "[Na+].[Na+].[O-][Se]([O-])(=O)=O",
      "sodium stibogluconate": "O.O.O.O.O.O.O.O.O.[Na+].[Na+].[Na+].[H][C@@]1(O[Sb]2(O)(O[C@@H](C([O-])=O)[C@@]1([H])O2)O[Sb]12([O-])O[C@@H](C([O-])=O)[C@@]([H])(O1)[C@]([H])(O2)[C@H](O)CO)[C@H](O)CO",
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")