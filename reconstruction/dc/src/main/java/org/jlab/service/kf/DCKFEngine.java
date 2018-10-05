package org.jlab.service.kf;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.track.TrackCandListFinder;

/**
 * An engine to run the Kalman filtering process for the DC software.
 *
 * @author benkel
 */
public class DCKFEngine extends ReconstructionEngine {

    DCGeant4Factory dcDetector;
    private int dataeventCounter;

    public DCKFEngine() {
        super("DCKF", "benkel", "0.1");
    }

    @Override
    public boolean init() {
        // super.LoadTables();
        Constants.Load();
        dataeventCounter = 0;
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        dataeventCounter++;
        System.out.println("[DCKF] Summary for DataEvent number " +
                           dataeventCounter + ":");

        // if (event.hasBank("MC::Particle") && this.getEngineConfigString("wireDistort")==null) {
        //     Constants.setWIREDIST(0);
        // }

        if (!event.hasBank("RUN::config")) return true;
        DataBank headerBank = event.getBank("RUN::config");

        if (headerBank.getInt("run", 0) == 0) return true;

        // TEST: GET org.jlab.rec.dc.Constants.HITBASE
        TrackCandListFinder trkCandFinder = new TrackCandListFinder(Constants.HITBASE);

        // TODO: GET org.jlab.rec.dc.cross.CrossList
        //       Maybe figure out a way to pass it through HBCrosses
        if (!event.hasBank("HitBasedTrkg::HBCrossLists"))
            System.out.println("        HBCrossLists cannot be found in the DataEvent.");
        else
            System.out.println("        HBCrossLists was found in the DataEvent!");

        // TEST: GET org.jlab.detector.geant4.v2.DCGeant4Factory
        // TEST: GET org.jlab.clas.swimtools.Swimmer.TORSCALE
        // TEST: GET org.jlab.clas.swimtools.Swim (magnetic field)
        Swim dcSwim = new Swim();


        // TODO: RUN trkCandFinder.getTrackCands()
        /* trkCandFinder.getTrackCands
        (
            ,
            dcDetector,
            Swimmer.getTorScale(),
            dcSwim
        );

        */

        // TODO: WRITE list<org.jlab.rec.dc.track.Track> to evio event

        System.out.println("        Processed DataEvent object " +
                           dataeventCounter + " successfully!\n");
        return true;
    }
}
