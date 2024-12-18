//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B1
{

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	PrimaryGeneratorAction::PrimaryGeneratorAction()
	{
		G4int n_particle = 1;
		fParticleGun = new G4ParticleGun(n_particle);

		// default particle kinematic
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
		G4String particleName;
		G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "gamma");
		fParticleGun->SetParticleDefinition(particle);
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
		fParticleGun->SetParticleEnergy(0.6 * MeV);
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	PrimaryGeneratorAction::~PrimaryGeneratorAction()
	{
		delete fParticleGun;
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
	{
		// Vị trí phát hạt (có thể giữ nguyên)
		G4double x0 = 0.0 * m;
		G4double y0 = 0.25 * m;
		G4double z0 = 0.0 * m;

		fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));

		// Tham số góc mở của nón (tính bằng radian)
		G4double theta_max = 20.0 * deg;  // Góc mở của nón, ví dụ 20 độ

		// Sinh góc phương vị (theta) trong khoảng [0, 2π]
		G4double theta = 2 * M_PI * G4UniformRand();  // Góc phương vị
		// Sinh góc cực (phi) trong khoảng [0, theta_max]
		G4double phi = theta_max * G4UniformRand();   // Góc cực nằm trong [0, theta_max]

		// Tính toán các thành phần của vectơ hướng
		G4double ux = std::sin(phi) * std::cos(theta);  // Thành phần x của vectơ hướng
		G4double uy = std::cos(phi); // Thành phần y của vectơ hướng (lúc nào cũng âm)
		G4double uz = -std::sin(phi) * std::sin(theta);                    // Thành phần z của vectơ hướng

		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

		// Đặt năng lượng của hạt
		fParticleGun->SetParticleEnergy(0.6 * MeV);

		// Sinh ra sự kiện (phát hạt)
		fParticleGun->GeneratePrimaryVertex(event);
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
