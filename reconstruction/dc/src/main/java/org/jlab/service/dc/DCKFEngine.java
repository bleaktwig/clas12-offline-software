package org.jlab.service.dc;

/**
 * A CLARA engine to run the Kalman filtering process for the DC software.
 *
 * @author benkel
 */
public class DCKFEngine extends DCEngine {

    public DCKFEngine() {
        super("DCKF");
    }

    @Override
    public boolean init() {
        // super.LoadTables();
        System.out.println("DCKFEngine initiated correctly!");
        return true;
    }
}
