# 一：找能量最小帧

```
#!perl
#Layer P:poly
#Layer C:crystal
#Layer H:water
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

###################
# define
my $cinit=80;
###################


my $doc = $Documents{"SDS_Layer_lie.xtd"};
my $newStudyTable = Documents->New("InteractionEnergy.std");
my $calcSheet = $newStudyTable->ActiveSheet;

$calcSheet->ColumnHeading(0) = "Frame";
$calcSheet->ColumnHeading(1) = "Cell";
$calcSheet->ColumnHeading(2) = "Total Energy of Cell";
$calcSheet->ColumnHeading(3) = "Layer P";
$calcSheet->ColumnHeading(4) = "Energy of LayerP";
$calcSheet->ColumnHeading(5) = "Layer C"; 
$calcSheet->ColumnHeading(6) = "Energy of LayerC";
$calcSheet->ColumnHeading(7) = "Interaction Energy";
$calcSheet->ColumnHeading(8) = "Interaction Energy per angstrom^2";


my $forcite = Modules->Forcite;
$forcite->ChangeSettings(Settings(CurrentForcefield => "COMPASS", WriteLevel => "Silent"));
my $lengthA = $doc->Lattice3D->LengthA;
my $lengthB = $doc->Lattice3D->LengthB;
my $surfaceArea = $lengthA * $lengthB;
print "The surface area is $surfaceArea angstrom^2\n";

my $numFrames = $doc->Trajectory->NumFrames;
print "Number of frames being analyzed is $numFrames\n";


my $dc=$cinit;

for (my $counter = $cinit ; $counter <= $numFrames; ++$counter) {
$doc->Trajectory->CurrentFrame = $counter;
my $atoms=$doc->UnitCell->Sets("Layer C")->Atoms;
$atoms->Unfix("XYZ");
my $allDoc = Documents->New("all.xsd");
my $layerpDoc = Documents->New("layerp.xsd");
my $layercDoc = Documents->New("layerc.xsd");
#relax carbon
$allDoc->CopyFrom($doc);
$allDoc->UnitCell->Sets("Layer H")->Atoms->Delete;
$layerpDoc->CopyFrom($doc);
$layerpDoc->UnitCell->Sets("Layer C")->Atoms->Delete;
$layerpDoc->UnitCell->Sets("Layer H")->Atoms->Delete;
$layercDoc->CopyFrom($doc);
$layercDoc->UnitCell->Sets("Layer P")->Atoms->Delete;
$layercDoc->UnitCell->Sets("Layer H")->Atoms->Delete;
$calcSheet->Cell($counter-$dc,0) = "No."+$counter;
$calcSheet->Cell($counter-$dc,1) = $allDoc;
$calcSheet->Cell($counter-$dc,3) = $layerpDoc;
$calcSheet->Cell($counter-$dc,5) = $layercDoc;
$forcite->Energy->Run($allDoc);
$calcSheet->Cell($counter-$dc, 2) = $allDoc->PotentialEnergy;
my $totalEnergy = $calcSheet->Cell($counter-$dc, 2);
$forcite->Energy->Run($layerpDoc);
$calcSheet->Cell($counter-$dc, 4) = $layerpDoc->PotentialEnergy;
my $layerpEnergy = $calcSheet->Cell($counter-$dc, 4);
$forcite->Energy->Run($layercDoc);
$calcSheet->Cell($counter-$dc, 6) = $layercDoc->PotentialEnergy;
my $layercEnergy = $calcSheet->Cell($counter-$dc, 6);
my $interactionEnergy = $totalEnergy - ($layerpEnergy + $layercEnergy);
$calcSheet->Cell($counter-$dc, 7) = $interactionEnergy;
my $interactionEnergyArea = $interactionEnergy / $surfaceArea;
$calcSheet->Cell($counter-$dc, 8) = $interactionEnergyArea;
$allDoc->Discard; 
$layerpDoc->Discard; 
$layercDoc->Discard; 
} 
```

# 二、氢键键长计算
```
#!perl

# Purpose: Calculate statistics (Min/Max/Mean) on the Hydrogen Bonds
#            in a structure. Put the results in a Study Table.

use strict;
use warnings;
use MaterialsScript qw(:all);

#Initialize variables for the stats calculations
my $totalLength = 0;
my $minLength   = 99999.9; #Arbitrary Big Num
my $maxLength   = 0;
my $row = 0;

#Get all the HBonds in the UnitCell
my $hbonds = $Documents{"D_T.xsd"}->UnitCell->HydrogenBonds;

#Create a new Study Table for the results
my $statsDoc = Documents->New("HBondStats.std");

foreach my $hbond (@$hbonds) {
    #Output the bond length for each HBond
    $statsDoc->Cell($row, 0) = "HBond $row";
    $statsDoc->Cell($row, 1) = $hbond->Length;

    #Update the statistics information
    $totalLength += $hbond->Length;

    if($hbond->Length < $minLength) {
        $minLength = $hbond->Length;
    }
    elsif($hbond->Length > $maxLength) {
        $maxLength = $hbond->Length;
    }

    ++$row;
}

#printout the overall statistics
$statsDoc->Cell($row, 0) = "Average";
$statsDoc->Cell($row, 1) = $totalLength/$row;
$statsDoc->Cell($row, 2) = "Min";
$statsDoc->Cell($row, 3) = $minLength;
$statsDoc->Cell($row, 4) = "Max";
$statsDoc->Cell($row, 5) = $maxLength;
```


# 三、原子键长计算
```
#!perl
#Layer P:poly
#Layer C:crystal
#Layer H:water
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

###################
# define
my $cinit=1;
my $step=10;
###################


my $doc = $Documents{"Gra_Layer_water_deleH.xtd"};
my $newStudyTable = Documents->New("HBondLength.std");
my $calcSheet = $newStudyTable->ActiveSheet;

$calcSheet->ColumnHeading(0) = "Frame";
$calcSheet->ColumnHeading(1) = "H1-O1";
$calcSheet->ColumnHeading(2) = "H2-O2";
$calcSheet->ColumnHeading(3) = "H3-O3";

my $dc=$cinit;
my $i=-1;

my $numFrames = $doc->Trajectory->NumFrames;
print "Number of frames being analyzed is $numFrames\n";

for (my $counter = $cinit; $counter <= $numFrames; $counter+=$step) {
$i+=1;
$doc->Trajectory->CurrentFrame = $counter;
my $myDoc = Documents->New("temp.xsd");
$myDoc->CopyFrom($doc);
$myDoc->NonPeriodicSuperstructure;
$calcSheet->Cell($counter-$dc-$i*$step+$i,0) = "Frame $counter";
my $atom=$myDoc->UnitCell->Sets("H1")->Atoms(0);
#my $atom=$myDoc->Atoms('O');
my $h_x=$atom->X;
my $h_y=$atom->Y;
my $h_z=$atom->Z;
#print "$h_x,\t$h_y,\t$h_z\n";
my $atom=$myDoc->UnitCell->Sets("O1")->Atoms(0);
my $O_x=$atom->X;
my $O_y=$atom->Y;
my $O_z=$atom->Z;
#print "$O_x,\t$O_y,\t$O_z\n";
my $l=sqrt(($h_x-$O_x)*($h_x-$O_x)+($h_y-$O_y)*($h_y-$O_y)+($h_z-$O_z)*($h_z-$O_z));
$calcSheet->Cell($counter-$dc-$i*$step+$i,1) = $l;
#print sqrt(2^2);

my $atom=$myDoc->UnitCell->Sets("H2")->Atoms(0);
#my $atom=$myDoc->Atoms('O');
my $h_x=$atom->X;
my $h_y=$atom->Y;
my $h_z=$atom->Z;
#print "$h_x,\t$h_y,\t$h_z\n";
my $atom=$myDoc->UnitCell->Sets("O2")->Atoms(0);
my $O_x=$atom->X;
my $O_y=$atom->Y;
my $O_z=$atom->Z;
#print "$h_x,\t$h_y,\t$h_z\n";
#print "$O_x,\t$O_y,\t$O_z\n";
my $l=sqrt(($h_x-$O_x)*($h_x-$O_x)+($h_y-$O_y)*($h_y-$O_y)+($h_z-$O_z)*($h_z-$O_z));
#print "$l\n";
$calcSheet->Cell($counter-$dc-$i*$step+$i,2) = $l;
my $atom=$myDoc->UnitCell->Sets("H3")->Atoms(0);
#my $atom=$myDoc->Atoms('O');
my $h_x=$atom->X;
my $h_y=$atom->Y;
my $h_z=$atom->Z;
#print "$h_x,\t$h_y,\t$h_z\n";
my $atom=$myDoc->UnitCell->Sets("O3")->Atoms(0);
my $O_x=$atom->X;
my $O_y=$atom->Y;
my $O_z=$atom->Z;
#print "$O_x,\t$O_y,\t$O_z\n";
my $l=sqrt(($h_x-$O_x)*($h_x-$O_x)+($h_y-$O_y)*($h_y-$O_y)+($h_z-$O_z)*($h_z-$O_z));
$calcSheet->Cell($counter-$dc-$i*$step+$i,3) = $l;
$myDoc->Discard; 

} 
```

# 四、Forcite实例
```
#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc = $Documents{"Layer.xtd"};
my $newStudyTable = Documents->New("InteractionEnergy.std");
my $calcSheet = $newStudyTable->ActiveSheet;

$calcSheet->ColumnHeading(0) = "Cell";
$calcSheet->ColumnHeading(1) = "Total Energy of Cell";
$calcSheet->ColumnHeading(2) = "Layer1";
$calcSheet->ColumnHeading(3) = "Energy of Layer1";
$calcSheet->ColumnHeading(4) = "Layer2"; 
$calcSheet->ColumnHeading(5) = "Energy of Layer2";
$calcSheet->ColumnHeading(6) = "Interaction Energy";
$calcSheet->ColumnHeading(7) = "Interaction Energy per angstrom^2";


my $forcite = Modules->Forcite;
$forcite->ChangeSettings(Settings(CurrentForcefield => "COMPASS", WriteLevel => "Silent"));
my $lengthA = $doc->Lattice3D->LengthA;
my $lengthB = $doc->Lattice3D->LengthB;
my $surfaceArea = $lengthA * $lengthB;
print "The surface area is $surfaceArea angstrom^2\n";

my $numFrames = $doc->Trajectory->NumFrames;
print "Number of frames being analyzed is $numFrames\n";

for (my $counter = 1; $counter <= 3; ++$counter) {
$doc->Trajectory->CurrentFrame = $counter;
my $allDoc = Documents->New("all.xsd");
my $layer1Doc = Documents->New("layer1.xsd");
my $layer2Doc = Documents->New("layer2.xsd");
$allDoc->CopyFrom($doc);
$layer1Doc->CopyFrom($doc);
$layer1Doc->UnitCell->Sets("Layer 2")->Atoms->Delete;
$layer2Doc->CopyFrom($doc);
$layer2Doc->UnitCell->Sets("Layer 1")->Atoms->Delete;
$calcSheet->Cell($counter-1,0) = $allDoc;
$calcSheet->Cell($counter-1,2) = $layer1Doc;
$calcSheet->Cell($counter-1,4) = $layer2Doc;
$forcite->Energy->Run($allDoc);
$calcSheet->Cell($counter-1, 1) = $allDoc->PotentialEnergy;
my $totalEnergy = $calcSheet->Cell($counter-1, 1);
$forcite->Energy->Run($layer1Doc);
$calcSheet->Cell($counter-1, 3) = $layer1Doc->PotentialEnergy;
my $layer1Energy = $calcSheet->Cell($counter-1, 3);
$forcite->Energy->Run($layer2Doc);
$calcSheet->Cell($counter-1, 5) = $layer2Doc->PotentialEnergy;
my $layer2Energy = $calcSheet->Cell($counter-1, 5);
my $interactionEnergy = $totalEnergy - ($layer1Energy + $layer2Energy);
$calcSheet->Cell($counter-1, 6) = $interactionEnergy;
my $interactionEnergyArea = $interactionEnergy / $surfaceArea;
$calcSheet->Cell($counter-1, 7) = $interactionEnergyArea;
$allDoc->Discard; 
$layer1Doc->Discard; 
$layer2Doc->Discard; 
} 
```
