package org.jlab.service.kf;

import org.jlab.io.base.DataEvent;
import org.jlab.clas.reco.ReconstructionEngine;
// import org.jlab.service.kf.KFTrackCandDataEvent;

/**
 * A CLARA engine to run the Kalman filtering process for the DC software.
 *
 * @author benkel
 */
public class DCKFEngine extends ReconstructionEngine {

    private int counter;

    public DCKFEngine() {
        super("DCKF", "benkel", "0.1");
    }

    @Override
    public boolean init() {
        // super.LoadTables();
        System.out.println("[DCKF] Initiated correctly!");
        counter = 0;
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        counter++;
        System.out.println("[DCKF] Processed DataEvent object " + counter + "!");
        return true;
    }

    public boolean processKFDataEvent(KFTrackCandDataEvent event) {
        System.out.println("[DCKF] Processed a KFTrackCandDataEvent object!");
        return true;
    }

    public boolean hello() {
        System.out.println("[DCKF] Am I even running?");
        return true;
    }
}
