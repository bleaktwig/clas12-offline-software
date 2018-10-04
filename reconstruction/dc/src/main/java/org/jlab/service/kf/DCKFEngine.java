package org.jlab.service.kf;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.reco.ReconstructionEngine;

import org.jlab.rec.dc.Constants;

/**
 * A CLARA engine to run the Kalman filtering process for the DC software.
 *
 * @author benkel
 */
public class DCKFEngine extends ReconstructionEngine {

    private int dataeventCounter;

    public DCKFEngine() {
        super("DCKF", "benkel", "0.1");
    }

    @Override
    public boolean init() {
        // super.LoadTables();
        System.out.println("[DCKF] Initiated correctly!");
        dataeventCounter = 0;
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        dataeventCounter++;
        System.out.println("[DCKF] Summary for DataEvent number " +
                           dataeventCounter + ":");

        Constants.Load();

        // if (event.hasBank("MC::Particle") && this.getEngineConfigString("wireDistort")==null) {
        //     Constants.setWIREDIST(0);
        // }

        if (!event.hasBank("RUN::config")) {
            System.out.println("[DCKF]   Header Bank is not in the DataEvent.");
            return true;
        }
        DataBank headerBank = event.getBank("RUN::config");

        if (headerBank.getInt("run", 0) == 0) {
            System.out.println("[DCKF]   Header Bank's run is equal to cero.");
            return true;
        }

        // GET org.jlab.rec.dc.Constants.HITBASE
        TrackCandListFinder trkcandFinder = new TrackCandListFinder(Constants.HITBASE);

        // TODO: GET org.jlab.rec.dc.cross.CrossList
        if (!event.hasBank("HitBasedTrkg::HBCrossLists")) {
            System.out.println("[DCKF]   HBCrossLists cannot be found in the DataEvent.");
        }
        else {
            System.out.println("[DCKF]   HBCrossLists was found in the DataEvent!");
        }

        // TODO: GET org.jlab.detector.geant4.v2.DCGeant4Factory
        // TODO: GET org.jlab.clas.swimtools.Swimmer.TORSCALE
        // TODO: GET org.jlab.clas.swimtools.Swim (instance)

        // TODO: WRITE list<org.jlab.rec.dc.track.Track>

        System.out.println("[DCKF]   Processed DataEvent object " +
                           dataeventCounter + " successfully!\n");
        return true;
    }
}
