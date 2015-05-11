/*	{
		DynamicList<label> overFilledCells;
		forAll(alpha1,ci)
		{
			if (alpha1[ci] > 1)
			{
				overFilledCells.append(ci);
			}
		}
		overFilledCells.shrink();
		label nOverFilledCells = overFilledCells.size();		
		label ci = -1;
		while (ci < nOverFilledCells-1)
		{
			ci++;
			label cLabel = overFilledCells[ci];
			Info << "Cell " << cLabel << ":" << endl;
			labelList faceLabels = mesh.cells()[cLabel];
			scalar magdV = 0;
			DynamicList<label> outFaces;
			forAll(faceLabels,fi)
			{
				label fLabel = faceLabels[fi];
				if (fLabel<mesh.nInternalFaces()) //face is internal and hence between two cells
				{
					if (cLabel == mesh.owner()[fLabel]) //positive flux means out of cell cLabel
					{
						if (phi[fLabel] > 0 &&  mag(dVf[fLabel]) < mag(phi[fLabel])*runTime.deltaT().value() )
						{
							magdV += mag(dVf[fLabel]);
							outFaces.append(fLabel);
						}
					}
					else //cLabel must be neighbour and negative phi means positive outflow 
					{
						if (phi[fLabel] < 0 && mag(dVf[fLabel]) < mag(phi[fLabel])*runTime.deltaT().value() )
						{
							magdV += mag(dVf[fLabel]);
							outFaces.append(fLabel);
						}
					}
				}
			}
			outFaces.shrink();
			if (magdV > SMALL)
			{
				forAll(outFaces,fi)
				{
					label fLabel = outFaces[fi];
					label otherCellLabel;
					if (cLabel == mesh.owner()[fLabel])
					{
						otherCellLabel = mesh.neighbour()[fLabel];
					}
					else
					{
						otherCellLabel = mesh.owner()[fLabel];					
					}
					scalar alphaOtherOld = alpha1[otherCellLabel];
					scalar dVi = (alpha1[cLabel]-1)*mesh.V()[cLabel]*mag(dVf[fLabel])/magdV;
					Info << "Transferring " << dVi << " m3 of water from cell " << cLabel << " to cell " << otherCellLabel << endl;
					Info << "Old alpha1[cLabel]: " <<  alpha1[cLabel] << ", old alpha[otherCellLabel] = " << alpha1[otherCellLabel] << endl;
					alpha1[otherCellLabel] += dVi/mesh.V()[otherCellLabel];
					alpha1[cLabel] -= dVi/mesh.V()[cLabel];
					Info << "New alpha1[cLabel]: " <<  alpha1[cLabel] << ", new alpha[otherCellLabel] = " << alpha1[otherCellLabel] << endl;
					if (alphaOtherOld<1 && alpha1[otherCellLabel]>1)
					{
						overFilledCells.append(otherCellLabel);
						nOverFilledCells++;
						Info << "Cell " << otherCellLabel << " got overfilled!!!" << endl;
					}
				}
			}
		}
	}
*/	
/////////////////////////////////////////

	Info << "-------------------------------" << endl;
	{
		DynamicList<label> overEmptiedCells;
		forAll(alpha1,ci)
		{
			if (alpha1[ci] < 0)
			{
				overEmptiedCells.append(ci);
			}
		}
		overEmptiedCells.shrink();
		label nOverEmptiedCells = overEmptiedCells.size();
		label ci = -1;
		while (ci < nOverEmptiedCells-1)
		{
			ci++;
			label cLabel = overEmptiedCells[ci];
			Info << "Cell " << cLabel << ":" << endl;
			labelList faceLabels = mesh.cells()[cLabel];
			scalar magdV = 0;
			DynamicList<label> outFaces;
			forAll(faceLabels,fi)
			{
				label fLabel = faceLabels[fi];
				if (fLabel<mesh.nInternalFaces()) //face is internal and hence between two cells
				{
					if (cLabel == mesh.owner()[fLabel]) //positive flux means out of cell cLabel
					{
						if (phi[fLabel] > 0)
						{
							magdV += mag(dVf[fLabel]);
							outFaces.append(fLabel);
						}
					}
					else
					{
						if (phi[fLabel] < 0)
						{
							magdV += mag(dVf[fLabel]);
							outFaces.append(fLabel);
						}
					}
				}
			}
			outFaces.shrink();
			
			if(magdV > SMALL)
			{
				forAll(outFaces,fi)
				{
					label fLabel = outFaces[fi];
					label otherCellLabel;
					if (cLabel == mesh.owner()[fLabel])
					{
						otherCellLabel = mesh.neighbour()[fLabel];
					}
					else
					{
						otherCellLabel = mesh.owner()[fLabel];					
					}
					scalar alphaOtherOld = alpha1[otherCellLabel];
					scalar dVi = (-alpha1[cLabel])*mesh.V()[cLabel]*mag(dVf[fLabel])/magdV;
					Info << "Transferring " << dVi << " m3 of water from cell " << otherCellLabel << " to cell " << cLabel << endl;
					Info << "Old alpha1[cLabel]: " <<  alpha1[cLabel] << ", old alpha[otherCellLabel] = " << alpha1[otherCellLabel] << endl;
					alpha1[otherCellLabel] -= dVi/mesh.V()[otherCellLabel];
					alpha1[cLabel] += dVi/mesh.V()[cLabel];
					Info << "New alpha1[cLabel]: " <<  alpha1[cLabel] << ", new alpha[otherCellLabel] = " << alpha1[otherCellLabel] << endl;
					if (alphaOtherOld>0 && alpha1[otherCellLabel]<0)
					{
						overEmptiedCells.append(otherCellLabel);
						nOverEmptiedCells++;
						Info << "Cell " << otherCellLabel << " got overemptied!!!" << endl;
					}
				}
			}
		}
	}