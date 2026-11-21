-- Generated from the pinned RDKit QED definitions by generate_fragments.py.
-- See RDKIT-DATA-LICENSE.txt.
import Descriptors.Fragment

namespace Descriptors.Fragment

-- [o&H0&X2]
private def qedAcceptors0 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 8 && (g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (2))), []⟩]

-- [O&H1&X2&v2]
private def qedAcceptors1 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))) && ((g.atom v).totalDegree == (2))) && ((g.atom v).valence == (2))), []⟩]

-- [O&H0&X2&v2]
private def qedAcceptors2 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (2))) && ((g.atom v).valence == (2))), []⟩]

-- [O&H0&X1&v2]
private def qedAcceptors3 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (1))) && ((g.atom v).valence == (2))), []⟩]

-- [O&-&X1]
private def qedAcceptors4 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))) && ((g.atom v).totalDegree == (1))), []⟩]

-- [S&H0&X2&v2]
private def qedAcceptors5 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (2))) && ((g.atom v).valence == (2))), []⟩]

-- [S&H0&X1&v2]
private def qedAcceptors6 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (1))) && ((g.atom v).valence == (2))), []⟩]

-- [S&-&X1]
private def qedAcceptors7 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))) && ((g.atom v).totalDegree == (1))), []⟩]

-- [n&H0&X2]
private def qedAcceptors8 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (2))), []⟩]

-- [N&H0&X1&v3]
private def qedAcceptors9 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (0))) && ((g.atom v).totalDegree == (1))) && ((g.atom v).valence == (3))), []⟩]

-- [$([N&+0&X3&v3])&!$(N[C,S]=O)]
private def qedAcceptors10 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => (((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).valence == (3))), []⟩]
  let sub1 : List Step := [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]
  [
    ⟨fun v => ((rooted g v sub0) && !(rooted g v sub1)), []⟩]

def qedAcceptors (g : Graph) : List (List Step) :=
  [qedAcceptors0 g,
   qedAcceptors1 g,
   qedAcceptors2 g,
   qedAcceptors3 g,
   qedAcceptors4 g,
   qedAcceptors5 g,
   qedAcceptors6 g,
   qedAcceptors7 g,
   qedAcceptors8 g,
   qedAcceptors9 g,
   qedAcceptors10 g]

-- *1[O,S,N]*1
private def qedAlerts0 (g : Graph) : List Step :=
  [
    ⟨fun _ => (true), []⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun _ => (true), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [S,C](=[O,S])[F,Br,Cl,I]
private def qedAlerts1 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) || ((g.atom v).z == 6 && !(g.atom v).aromatic)), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&X4][Cl,Br,I]
private def qedAlerts2 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), []⟩,
    ⟨fun v => ((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#6]S(=O)(=O)O[#6]
private def qedAlerts3 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [$([C&H1]),$(CC)]#CC(=O)[#6]
private def qedAlerts4 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩]
  let sub1 : List Step := [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  [
    ⟨fun v => ((rooted g v sub0) || (rooted g v sub1)), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [$([C&H1]),$(CC)]#CC(=O)O[#6]
private def qedAlerts5 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩]
  let sub1 : List Step := [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  [
    ⟨fun v => ((rooted g v sub0) || (rooted g v sub1)), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- n[O&H1]
private def qedAlerts6 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [$([C&H1]),$(CC)]#CS(=O)(=O)[#6]
private def qedAlerts7 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩]
  let sub1 : List Step := [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  [
    ⟨fun v => ((rooted g v sub0) || (rooted g v sub1)), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C=C(C=O)C=O
private def qedAlerts8 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- n1c([F,Cl,Br,I])cccc1
private def qedAlerts9 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&H1]=O
private def qedAlerts10 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [#8][#8]
private def qedAlerts11 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (8)), []⟩,
    ⟨fun v => ((g.atom v).z == (8)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&!R]=[N&!R]
private def qedAlerts12 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), []⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), [(0, fun b => (b % 32 == 2))]⟩]

-- [N&!R]=[N&!R]
private def qedAlerts13 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), []⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), [(0, fun b => (b % 32 == 2))]⟩]

-- [#6](=O)[#6]=O
private def qedAlerts14 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [#16][#16]
private def qedAlerts15 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (16)), []⟩,
    ⟨fun v => ((g.atom v).z == (16)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#7][N&H2]
private def qedAlerts16 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (7)), []⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C(=O)N[N&H2]
private def qedAlerts17 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#6]=S
private def qedAlerts18 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [$([C&H2]),$([C&H1][C&X4]),$(C([C&X4])[C&X4])]=[$([C&H2]),$([C&H1][C&X4]),$(C([C&X4])[C&X4])]
private def qedAlerts19 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), []⟩]
  let sub1 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  let sub2 : List Step := [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  let sub3 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), []⟩]
  let sub4 : List Step := [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  let sub5 : List Step := [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  [
    ⟨fun v => (((rooted g v sub0) || (rooted g v sub1)) || (rooted g v sub2)), []⟩,
    ⟨fun v => (((rooted g v sub3) || (rooted g v sub4)) || (rooted g v sub5)), [(0, fun b => (b % 32 == 2))]⟩]

-- C1(=[O,N])C=CC(=[O,N])C=C1
private def qedAlerts20 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C1(=[O,N])C(=[O,N])C=CC=C1
private def qedAlerts21 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- a12aa3a(aa1aaaa2)aaaa3
private def qedAlerts22 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (8, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (10, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- a12a(a3a(aa1)aaaa3)aaaa2
private def qedAlerts23 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(8, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (12, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- a1aa2a3a(a1)A=AA=A3=AA=A2
private def qedAlerts24 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 2)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (!(g.atom v).aromatic), [(0, fun b => (b % 32 == 2)), (9, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1cc([N&H2])ccc1
private def qedAlerts25 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#80,#26,#33,#51,#30,#34,se,#52,B,#14,#11,#20,#32,#47,#12,#19,#56,#38,#4,#22,#42,#25,#44,#46,#28,#29,#79,#48,#13,#31,#50,#45,#81,#83,#41,#3,#82,#72,#67]
private def qedAlerts26 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((((((((((((((((((((((((((((((((((((((g.atom v).z == (80)) || ((g.atom v).z == (26))) || ((g.atom v).z == (33))) || ((g.atom v).z == (51))) || ((g.atom v).z == (30))) || ((g.atom v).z == (34))) || ((g.atom v).z == 34 && (g.atom v).aromatic)) || ((g.atom v).z == (52))) || ((g.atom v).z == 5 && !(g.atom v).aromatic)) || ((g.atom v).z == (14))) || ((g.atom v).z == (11))) || ((g.atom v).z == (20))) || ((g.atom v).z == (32))) || ((g.atom v).z == (47))) || ((g.atom v).z == (12))) || ((g.atom v).z == (19))) || ((g.atom v).z == (56))) || ((g.atom v).z == (38))) || ((g.atom v).z == (4))) || ((g.atom v).z == (22))) || ((g.atom v).z == (42))) || ((g.atom v).z == (25))) || ((g.atom v).z == (44))) || ((g.atom v).z == (46))) || ((g.atom v).z == (28))) || ((g.atom v).z == (29))) || ((g.atom v).z == (79))) || ((g.atom v).z == (48))) || ((g.atom v).z == (13))) || ((g.atom v).z == (31))) || ((g.atom v).z == (50))) || ((g.atom v).z == (45))) || ((g.atom v).z == (81))) || ((g.atom v).z == (83))) || ((g.atom v).z == (41))) || ((g.atom v).z == (3))) || ((g.atom v).z == (82))) || ((g.atom v).z == (72))) || ((g.atom v).z == (67))), []⟩]

-- I
private def qedAlerts27 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 53 && !(g.atom v).aromatic), []⟩]

-- OS(=O)(=O)[O&-]
private def qedAlerts28 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [N&+](=O)[O&-]
private def qedAlerts29 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C(=O)N[O&H1]
private def qedAlerts30 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C1NC(=O)NC1=O
private def qedAlerts31 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [S&H1]
private def qedAlerts32 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩]

-- [S&-]
private def qedAlerts33 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))), []⟩]

-- c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]
private def qedAlerts34 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]
private def qedAlerts35 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 17 && !(g.atom v).aromatic) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 9 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&R1]1[C&R1][C&R1][C&R1][C&R1][C&R1][C&R1]1
private def qedAlerts36 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&R1]1[C&R1][C&R1]cc[C&R1][C&R1]1
private def qedAlerts37 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&R2]1[C&R2][C&R2][C&R2][C&R2][C&R2][C&R2][C&R2]1
private def qedAlerts38 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&R2]1[C&R2][C&R2]cc[C&R2][C&R2][C&R2]1
private def qedAlerts39 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&H2&R2]1N[C&H2&R2][C&H2&R2][C&H2&R2][C&H2&R2][C&H2&R2]1
private def qedAlerts40 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&H2&R2]1N[C&H2&R2][C&H2&R2][C&H2&R2][C&H2&R2][C&H2&R2][C&H2&R2]1
private def qedAlerts41 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C#C
private def qedAlerts42 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- [O&R2,N&R2]@[C&R2]@[C&R2]@[O&R2,N&R2]@[C&R2]@[C&R2]@[O&R2,N&R2]
private def qedAlerts43 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).rings == (2)))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b >= 32))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b >= 32))]⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).rings == (2)))), [(0, fun b => (b >= 32))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b >= 32))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b >= 32))]⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (2))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).rings == (2)))), [(0, fun b => (b >= 32))]⟩]

-- [$([N&+&R]),$([n&+&R]),$([N&+]=C)][O&-]
private def qedAlerts44 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))) && ((g.atom v).rings != 0)), []⟩]
  let sub1 : List Step := [
    ⟨fun v => ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).charge == (1))) && ((g.atom v).rings != 0)), []⟩]
  let sub2 : List Step := [
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]
  [
    ⟨fun v => (((rooted g v sub0) || (rooted g v sub1)) || (rooted g v sub2)), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#6]=N[O&H1]
private def qedAlerts45 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#6]=NOC=O
private def qedAlerts46 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [#6](=O)[C&X4,C&R0&X3,O][#6]=O
private def qedAlerts47 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))) || ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))) && ((g.atom v).totalDegree == (3)))) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- c1ccc2c(c1)ccc(=O)o2
private def qedAlerts48 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [O&+,o&+,S&+,s&+]
private def qedAlerts49 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))) || (((g.atom v).z == 8 && (g.atom v).aromatic) && ((g.atom v).charge == (1)))) || (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).charge == (1)))) || (((g.atom v).z == 16 && (g.atom v).aromatic) && ((g.atom v).charge == (1)))), []⟩]

-- N=C=O
private def qedAlerts50 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [N&X3,N&X4][F,Cl,Br,I]
private def qedAlerts51 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4)))), []⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1ccccc1OC(=O)[#6]
private def qedAlerts52 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&R0]=[C&R0][C&R0]=[C&R0]
private def qedAlerts53 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 2))]⟩]

-- [C&+,c&+,C&-,c&-]
private def qedAlerts54 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))) || (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).charge == (1)))) || (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1)))) || (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).charge == (-1)))), []⟩]

-- N=[N&+]=[N&-]
private def qedAlerts55 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))), [(0, fun b => (b % 32 == 2))]⟩]

-- C12C(NC(N1)=O)CSC2
private def qedAlerts56 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (3, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (7, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1c([O&H1])c([O&H1,N&H2,N&H1])ccc1
private def qedAlerts57 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2)))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- P
private def qedAlerts58 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 15 && !(g.atom v).aromatic), []⟩]

-- [N,O,S]C#N
private def qedAlerts59 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 7 && !(g.atom v).aromatic) || ((g.atom v).z == 8 && !(g.atom v).aromatic)) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- C=C=O
private def qedAlerts60 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- [#14][F,Cl,Br,I]
private def qedAlerts61 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (14)), []⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [S&X2]O
private def qedAlerts62 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#14&R0,C&R0](c1ccccc1)(c1ccccc1)c1ccccc1
private def qedAlerts63 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == (14)) && ((g.atom v).rings == (0))) || (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0)))), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(12, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- O1CCCCC1OC1CCC2CCCCC2C1
private def qedAlerts64 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (8, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- N=[C&R0][N,n,O,S]
private def qedAlerts65 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => (((((g.atom v).z == 7 && !(g.atom v).aromatic) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [c&R2]1[c&R2][c&R2]([N&v3&X3,N&v4&X4])[c&R2][c&R2][c&R2]1[c&R2]1[c&R2][c&R2][c&R2]([N&v3&X3,N&v4&X4])[c&R2][c&R2]1
private def qedAlerts66 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).valence == (3))) && ((g.atom v).totalDegree == (3))) || ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).valence == (4))) && ((g.atom v).totalDegree == (4)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).valence == (3))) && ((g.atom v).totalDegree == (3))) || ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).valence == (4))) && ((g.atom v).totalDegree == (4)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C=[C&!R]C#N
private def qedAlerts67 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings > 0)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- [c&R2]1[c&R2]c([N&+0&X3&R0,n&X3&R0])c([N&+0&X3&R0,n&X3&R0])[c&R2][c&R2]1
private def qedAlerts68 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [c&R2]1[c&R2]c([N&+0&X3&R0,n&X3&R0])[c&R2]c([N&+0&X3&R0,n&X3&R0])[c&R2]1
private def qedAlerts69 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (6, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [c&R2]1[c&R2]c([N&+0&X3&R0,n&X3&R0])[c&R2][c&R2]c1[N&+0&X3&R0,n&X3&R0]
private def qedAlerts70 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [O&H1]c1ccc([O&H1,N&H2,N&H1])cc1
private def qedAlerts71 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2)))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1ccccc1OC(=O)O
private def qedAlerts72 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [S&X2&H0]N
private def qedAlerts73 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))) && ((g.atom v).hydrogens == (0))), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c12ccccc1SC(S)=N2
private def qedAlerts74 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2)), (8, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c12ccccc1SC(=S)N2
private def qedAlerts75 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (8, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1nnnn1C=O
private def qedAlerts76 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (3, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- s1c(S)nnc1NC=O
private def qedAlerts77 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 16 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- S1C=CSC1=S
private def qedAlerts78 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (3, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- C(=O)Onnn
private def qedAlerts79 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- OS(=O)(=O)C(F)(F)F
private def qedAlerts80 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- N#CC[O&H1]
private def qedAlerts81 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- N#CC=O
private def qedAlerts82 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- S(=O)(=O)C#N
private def qedAlerts83 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- N[C&H2]C#N
private def qedAlerts84 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- C1(=O)NCC1
private def qedAlerts85 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (3, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- S(=O)(=O)[O&-,O&H1]
private def qedAlerts86 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))) || (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1)))), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- NC[F,Cl,Br,I]
private def qedAlerts87 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C=[C&!R]O
private def qedAlerts88 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings > 0)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [N&X2&+0]=[O&+0]
private def qedAlerts89 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))) && ((g.atom v).charge == (0))), []⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (0))), [(0, fun b => (b % 32 == 2))]⟩]

-- [O&R0,N&R0][O&R0,N&R0]
private def qedAlerts90 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).rings == (0)))), []⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]
private def qedAlerts91 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) || ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) || ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) || ((g.atom v).hydrogens == (1))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [C&X2&R0][N&X3&R0]
private def qedAlerts92 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))) && ((g.atom v).rings == (0))), []⟩,
    ⟨fun v => ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1ccccc1[C&!R]=[C&!R]c1ccccc1
private def qedAlerts93 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [N&X3&R0,N&X4&R0,O&R0,S&X2&R0][C&X4][N&X3&R0,N&X4&R0,O&R0,S&X2&R0]
private def qedAlerts94 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))) && ((g.atom v).rings == (0)))) || (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (0)))) || ((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))) && ((g.atom v).rings == (0)))), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).rings == (0))) || ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (4))) && ((g.atom v).rings == (0)))) || (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).rings == (0)))) || ((((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))) && ((g.atom v).rings == (0)))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [s,S,c,C,n,N,o,O]~[n&+,N&+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]
private def qedAlerts95 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)) || ((g.atom v).z == 8 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), []⟩,
    ⟨fun v => ((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).charge == (1))) || (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1)))), [(0, fun _ => (true))]⟩,
    ⟨fun v => (((((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)) || ((g.atom v).z == 8 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), [(0, fun _ => (true))]⟩,
    ⟨fun v => (((((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)) || ((g.atom v).z == 8 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), [(1, fun _ => (true))]⟩,
    ⟨fun v => (((((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)) || ((g.atom v).z == 8 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), [(2, fun _ => (true))]⟩]

-- [s,S,c,C,n,N,o,O]~[n&X3&+,N&X3&+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]
private def qedAlerts96 (g : Graph) : List Step :=
  [
    ⟨fun v => (((((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)) || ((g.atom v).z == 8 && (g.atom v).aromatic)) || ((g.atom v).z == 8 && !(g.atom v).aromatic)), []⟩,
    ⟨fun v => (((((g.atom v).z == 7 && (g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).charge == (1))) || ((((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))) && ((g.atom v).charge == (1)))), [(0, fun _ => (true))]⟩,
    ⟨fun v => (((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(0, fun _ => (true))]⟩,
    ⟨fun v => (((((((g.atom v).z == 16 && (g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)) || ((g.atom v).z == 6 && (g.atom v).aromatic)) || ((g.atom v).z == 6 && !(g.atom v).aromatic)) || ((g.atom v).z == 7 && (g.atom v).aromatic)) || ((g.atom v).z == 7 && !(g.atom v).aromatic)), [(1, fun _ => (true))]⟩]

-- *=[N&+]=*
private def qedAlerts97 (g : Graph) : List Step :=
  [
    ⟨fun _ => (true), []⟩,
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).charge == (1))), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun _ => (true), [(0, fun b => (b % 32 == 2))]⟩]

-- [S&X3](=O)[O&-,O&H1]
private def qedAlerts98 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (3))), []⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).charge == (-1))) || (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).hydrogens == (1)))), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- N#N
private def qedAlerts99 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 3))]⟩]

-- F.F.F.F
private def qedAlerts100 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), []⟩]

-- [R0&D2][R0&D2][R0&D2][R0&D2]
private def qedAlerts101 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).rings == (0)) && ((g.atom v).degree == (2))), []⟩,
    ⟨fun v => (((g.atom v).rings == (0)) && ((g.atom v).degree == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).rings == (0)) && ((g.atom v).degree == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).rings == (0)) && ((g.atom v).degree == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [c&R,C&R]~C(=O)NC(=O)~[c&R,C&R]
private def qedAlerts102 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings != 0)) || (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings != 0))), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun _ => (true))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((((g.atom v).z == 6 && (g.atom v).aromatic) && ((g.atom v).rings != 0)) || (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings != 0))), [(1, fun _ => (true))]⟩]

-- C=&!@CC=[O,S]
private def qedAlerts103 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => ((b % 32 == 2) && !(b >= 32)))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩]

-- [#6,#8,#16][#6](=O)O[#6]
private def qedAlerts104 (g : Graph) : List Step :=
  [
    ⟨fun v => ((((g.atom v).z == (6)) || ((g.atom v).z == (8))) || ((g.atom v).z == (16))), []⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c[C&R0](=[O,S])[#6]
private def qedAlerts105 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).rings == (0))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == (6)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c[S&X2][C&!R]
private def qedAlerts106 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).totalDegree == (2))), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && !((g.atom v).rings != 0)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- C=C=C
private def qedAlerts107 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩]

-- c1nc([F,Cl,Br,I,S])ncc1
private def qedAlerts108 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1ncnc([F,Cl,Br,I,S])c1
private def qedAlerts109 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)) || ((g.atom v).z == 16 && !(g.atom v).aromatic)), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12)), (5, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- c1ncc2c(n1)nc(n2)[F,Cl,Br,I]
private def qedAlerts110 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), []⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 7 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => (((((g.atom v).z == 9 && !(g.atom v).aromatic) || ((g.atom v).z == 17 && !(g.atom v).aromatic)) || ((g.atom v).z == 35 && !(g.atom v).aromatic)) || ((g.atom v).z == 53 && !(g.atom v).aromatic)), [(1, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [#6]S(=O)(=O)c1ccc(cc1)F
private def qedAlerts111 (g : Graph) : List Step :=
  [
    ⟨fun v => ((g.atom v).z == (6)), []⟩,
    ⟨fun v => ((g.atom v).z == 16 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(0, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 8 && !(g.atom v).aromatic), [(1, fun b => (b % 32 == 2))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 6 && (g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12)), (4, fun b => (b % 32 == 1 || b % 32 == 12))]⟩,
    ⟨fun v => ((g.atom v).z == 9 && !(g.atom v).aromatic), [(2, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]

-- [N&15*]
private def qedAlerts112 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 7 && !(g.atom v).aromatic) && ((g.atom v).isotope == (15))), []⟩]

-- [C&13*]
private def qedAlerts113 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 6 && !(g.atom v).aromatic) && ((g.atom v).isotope == (13))), []⟩]

-- [O&18*]
private def qedAlerts114 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 8 && !(g.atom v).aromatic) && ((g.atom v).isotope == (18))), []⟩]

-- [S&34*]
private def qedAlerts115 (g : Graph) : List Step :=
  [
    ⟨fun v => (((g.atom v).z == 16 && !(g.atom v).aromatic) && ((g.atom v).isotope == (34))), []⟩]

def qedAlerts (g : Graph) : List (List Step) :=
  [qedAlerts0 g,
   qedAlerts1 g,
   qedAlerts2 g,
   qedAlerts3 g,
   qedAlerts4 g,
   qedAlerts5 g,
   qedAlerts6 g,
   qedAlerts7 g,
   qedAlerts8 g,
   qedAlerts9 g,
   qedAlerts10 g,
   qedAlerts11 g,
   qedAlerts12 g,
   qedAlerts13 g,
   qedAlerts14 g,
   qedAlerts15 g,
   qedAlerts16 g,
   qedAlerts17 g,
   qedAlerts18 g,
   qedAlerts19 g,
   qedAlerts20 g,
   qedAlerts21 g,
   qedAlerts22 g,
   qedAlerts23 g,
   qedAlerts24 g,
   qedAlerts25 g,
   qedAlerts26 g,
   qedAlerts27 g,
   qedAlerts28 g,
   qedAlerts29 g,
   qedAlerts30 g,
   qedAlerts31 g,
   qedAlerts32 g,
   qedAlerts33 g,
   qedAlerts34 g,
   qedAlerts35 g,
   qedAlerts36 g,
   qedAlerts37 g,
   qedAlerts38 g,
   qedAlerts39 g,
   qedAlerts40 g,
   qedAlerts41 g,
   qedAlerts42 g,
   qedAlerts43 g,
   qedAlerts44 g,
   qedAlerts45 g,
   qedAlerts46 g,
   qedAlerts47 g,
   qedAlerts48 g,
   qedAlerts49 g,
   qedAlerts50 g,
   qedAlerts51 g,
   qedAlerts52 g,
   qedAlerts53 g,
   qedAlerts54 g,
   qedAlerts55 g,
   qedAlerts56 g,
   qedAlerts57 g,
   qedAlerts58 g,
   qedAlerts59 g,
   qedAlerts60 g,
   qedAlerts61 g,
   qedAlerts62 g,
   qedAlerts63 g,
   qedAlerts64 g,
   qedAlerts65 g,
   qedAlerts66 g,
   qedAlerts67 g,
   qedAlerts68 g,
   qedAlerts69 g,
   qedAlerts70 g,
   qedAlerts71 g,
   qedAlerts72 g,
   qedAlerts73 g,
   qedAlerts74 g,
   qedAlerts75 g,
   qedAlerts76 g,
   qedAlerts77 g,
   qedAlerts78 g,
   qedAlerts79 g,
   qedAlerts80 g,
   qedAlerts81 g,
   qedAlerts82 g,
   qedAlerts83 g,
   qedAlerts84 g,
   qedAlerts85 g,
   qedAlerts86 g,
   qedAlerts87 g,
   qedAlerts88 g,
   qedAlerts89 g,
   qedAlerts90 g,
   qedAlerts91 g,
   qedAlerts92 g,
   qedAlerts93 g,
   qedAlerts94 g,
   qedAlerts95 g,
   qedAlerts96 g,
   qedAlerts97 g,
   qedAlerts98 g,
   qedAlerts99 g,
   qedAlerts100 g,
   qedAlerts101 g,
   qedAlerts102 g,
   qedAlerts103 g,
   qedAlerts104 g,
   qedAlerts105 g,
   qedAlerts106 g,
   qedAlerts107 g,
   qedAlerts108 g,
   qedAlerts109 g,
   qedAlerts110 g,
   qedAlerts111 g,
   qedAlerts112 g,
   qedAlerts113 g,
   qedAlerts114 g,
   qedAlerts115 g]

-- [$([A&R][!a])]
private def qedAliphatic0 (g : Graph) : List Step :=
  let sub0 : List Step := [
    ⟨fun v => ((!(g.atom v).aromatic) && ((g.atom v).rings != 0)), []⟩,
    ⟨fun v => !((g.atom v).aromatic), [(0, fun b => (b % 32 == 1 || b % 32 == 12))]⟩]
  [
    ⟨fun v => (rooted g v sub0), []⟩]

def qedAliphatic (g : Graph) : List (List Step) :=
  [qedAliphatic0 g]

end Descriptors.Fragment
