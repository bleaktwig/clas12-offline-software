package org.jlab.rec.dc.trajectory;

import java.util.ArrayList;
import java.util.List;

import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.detector.geant4.v2.FTOFGeant4Factory;
import org.jlab.detector.geant4.v2.ECGeant4Factory;
import org.jlab.detector.geant4.v2.PCALGeant4Factory;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;

/**
 * A class to load the geometry constants used in the DC reconstruction. The coordinate system used
 * in the tilted sector coordinate system.
 * @author ziegler
 */
public class TrajectorySurfaces {

    private List<ArrayList<Surface>> _DetectorPlanes = new ArrayList<ArrayList<Surface>>();

    double FVT_Z1stlayer  = 30.2967; // z-distance between target center and strips of the first
                                     //     layer.
    double FVT_Interlayer = 1.190;   // Keep this for now until the Geometry service is ready... or
                                     //     remove FMT from traj.

    public List<ArrayList<Surface>> getDetectorPlanes() {return _DetectorPlanes;}
    public synchronized void setDetectorPlanes(List<ArrayList<Surface>> aDetectorPlanes) {
        _DetectorPlanes = aDetectorPlanes;
    }

    public void LoadSurfaces(DCGeant4Factory dcDetector,
                             FTOFGeant4Factory ftofDetector,
                             ECGeant4Factory ecDetector,
                             PCALGeant4Factory pcalDetector) {

        int iw = 0;
        double d = 0;
        Vector3D n;

        for (int is = 0; is < 6; is++) {
            int index = 0;

            // System.out.println("CREATING SURFACES FOR SECTOR " + (is + 1));
            this._DetectorPlanes.add(new ArrayList<Surface>());

            // Add FMT
            for (int i = 0; i < 6; i++) {
                d = FVT_Z1stlayer + i*FVT_Interlayer;

                this._DetectorPlanes.get(is).add(new Surface("FMT" + (index+1), index++,
                                                             i+1, d, 0., 0., 1.));
            }
            index = 7; // End of MM + HTCC

            // Add DC
            n = this.RotateFromTSCtoLabC(0, 0, 1, is+1).toVector3D();
            for (int isup = 0; isup < 6; isup++) {
                for (int il = 0; il < 6; il++) {
                    d = dcDetector.getWireMidpoint(is, isup, il, iw).z;
                    this._DetectorPlanes.get(is).add(new Surface("DC" + (index - 6), index++,
                                                                 is*6 + il + 1, d,
                                                                 n.x(), n.y(), n.z()));
                }
            }
            index = 7; // End of MM + HTCC);

            // Outer detectors
            // LTCC
            this._DetectorPlanes.get(is).add(new Surface("LTCC", index++, 1, 624.23,
                                                         n.x(), n.y(), n.z()));
            index = 7; // End of MM + HTCC);

            // FTOF 2
            Vector3D P = ftofDetector.getFrontalFace(is+1, 3).point().toVector3D();
            n = ftofDetector.getFrontalFace(is+1, 3).normal();
            d = P.dot(n);
            this._DetectorPlanes.get(is).add(new Surface("FTOF2", index++, 3, -d,
                                                         -n.x(), -n.y(), -n.z()));
            index = 7; // End of MM + HTCC);

            // FTOF 18
            P = ftofDetector.getFrontalFace(is + 1, 2).point().toVector3D();
            n = ftofDetector.getFrontalFace(is + 1, 2).normal();
            d = P.dot(n);
            this._DetectorPlanes.get(is).add(new Surface("FTOF1b", index++, 2, -d,
                                                         -n.x(), -n.y(), -n.z()));
            index = 7; // End of MM + HTCC);

            // FTOF 1A
            P = ftofDetector.getFrontalFace(is + 1, 1).point().toVector3D();
            n = ftofDetector.getFrontalFace(is + 1, 1).normal();
            d = P.dot(n);
            this._DetectorPlanes.get(is).add(new Surface("FTOF1a", index++, 1, -d,
                                                         -n.x(), -n.y(), -n.z()));
            index = 7; // End of MM + HTCC);

            // PCAL (3)
            P = pcalDetector.getFrontalFace(is + 1).point().toVector3D();
            n = pcalDetector.getFrontalFace(is + 1).normal();
            d = P.dot(n);
            this._DetectorPlanes.get(is).add(new Surface("PCAL", index++, 1, -d,
                                                         -n.x(), -n.y(), -n.z()));
            index = 7; // End of MM + HTCC);

            // ECin (3)
            P = ecDetector.getFrontalFace(is + 1).point().toVector3D();
            n = ecDetector.getFrontalFace(is + 1).normal();
            d = P.dot(n);
            this._DetectorPlanes.get(is).add(new Surface("EC", index++, 1, -d,
                                                         -n.x(), -n.y(), -n.z()));
            index = 7; // End of MM + HTCC);
        }
    }

    private Point3D RotateFromTSCtoLabC(double X, double Y, double Z, int sector) {

        // Math.sin(Math.toRadians(25.)) = 0.42261826174069944
        // Math.cos(Math.toRadians(25.)) = 0.9063077870366499
        double rzs = -X * 0.42261826174069944 + Z * 0.9063077870366499;
        double rxs =  X * 0.9063077870366499  + Z * 0.42261826174069944;

        // Math.toRadians(60.) = 1.0471975511965976
        double rx = rxs * Math.cos((sector - 1) * 1.0471975511965976)
                    - Y * Math.sin((sector - 1) * 1.0471975511965976);
        double ry = rxs * Math.sin((sector - 1) * 1.0471975511965976)
                    + Y * Math.cos((sector - 1) * 1.0471975511965976);

        return new Point3D(rx, ry, rzs);
    }
}
