package org.jlab.rec.dc.hit;

import org.jlab.detector.geant4.v2.DCGeant4Factory;

/**
 * A DC hit characterized by superlayer, layer, sector, wire number, and time. The TDC to time
 * conversion has been done.
 *
 * @author ziegler
 */
public class Hit implements Comparable<Hit> {
    // Class implements Comparable interface to allow for sorting a collection of hits by wire
    //     number values

    private int _Id;
    private int _Sector;      // sector     [1...6]
    private int _Superlayer;  // superlayer [1...6]
    private int _Layer;    	  // layer      [1...6]
    private int _Wire;    	  // wire     [1...112]
    private int _TDC;         // NOTE: What does TDC stand for?
    private double _cellSize; // NOTE: What is this cell size?
    private double _DocaErr;  // Doca uncertainty, or error on the time in ns (4ns time window used
                              //     by default in reconstructing simulated data).

    public Hit(int sector, int superlayer, int layer, int wire, int TDC, int Id) {
        this._Sector = sector;
        this._Superlayer = superlayer;
        this._Layer = layer;
        this._Wire = wire;
        this._TDC = TDC;
        this._Id = Id;
    }

    public int get_Sector() {return _Sector;}
    public void set_Sector(int _Sector) {this._Sector = _Sector;}

    public int get_Superlayer() {return _Superlayer;}
    public void set_Superlayer(int _Superlayer) {this._Superlayer = _Superlayer;}

    public int get_Layer() {return _Layer;}
    public void set_Layer(int _Layer) {this._Layer = _Layer;}

    public int get_Wire() {return _Wire;}
    public void set_Wire(int _Wire) {this._Wire = _Wire;}

    public int get_TDC() {return _TDC;}
    public void set_TDC(int TDC) {this._TDC = TDC;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public int get_Region() {return (this._Superlayer + 1) / 2;}
    public int get_RegionSlayer() {return (this._Superlayer + 1) % 2 + 1;}

    public double get_CellSize() {return _cellSize;}
    public void set_CellSize(double cellSize) {_cellSize = cellSize;}

    // NOTE: Lacks JavaDoc comment
    public void calc_CellSize(DCGeant4Factory DcDetector) {
        double layerDiffAtMPln = DcDetector.getWireMidpoint(this.get_Sector() - 1,
                                                            this.get_Superlayer() - 1, 0, 0).x -
                                 DcDetector.getWireMidpoint(this.get_Sector() - 1,
                                                            this.get_Superlayer() - 1, 0, 1).x;

        _cellSize = 0.5 * Math.abs(layerDiffAtMPln);
    }

    public double get_DocaErr() {return _DocaErr;}
    public void set_DocaErr(double _docaErr) {this._DocaErr = _docaErr;}

    /**
     * Returns an int to sort a collection of hits by wire number. Sorting by wire is used in
     * clustering.
     * @param arg0 hit to compare to
     * @return     the sorting int
     */
    @Override
    public int compareTo (Hit arg) {
        int return_val = 0;
        int CompSec = this.get_Sector() < arg.get_Sector() ? -1 :
                      this.get_Sector() == arg.get_Sector() ? 0 : 1;
        int CompSly = this.get_Superlayer() < arg.get_Superlayer() ? -1 :
                      this.get_Superlayer() == arg.get_Superlayer() ? 0 : 1;
        int CompLay = this.get_Layer() < arg.get_Layer() ? -1 :
                      this.get_Layer() == arg.get_Layer() ? 0 : 1;
        int CompWir = this.get_Wire() < arg.get_Wire() ? -1 :
                      this.get_Wire() == arg.get_Wire() ? 0 : 1;

        int return_val1 = ((CompLay == 0) ? CompWir : CompLay);
        int return_val2 = ((CompSly == 0) ? return_val1 : CompSly);
        return ((CompSec == 0) ? return_val2 : CompSec);
    }

    /**
     * Calculates the center of the cell as a function of wire number in the local superlayer
     * coordinate system.
     * @param layer layer number from 1 to 6
     * @param wire wire number from 1 to 112
     */
    public double calcLocY(int layer, int wire) {
        // In old mc, layer 1 is closer to the beam than layer 2, while in hardware it's the
        //     opposite.
        // double brickwallPattern = GeometryLoader.getDcDetector().getWireMidpoint(0, 1, 1).x -
        //                           GeometryLoader.getDcDetector().getWireMidpoint(0, 0, 1).x;
        // double brickwallSign = Math.signum(brickwallPattern);
        double brickwallSign = -1;

        // center of the cell asfcn wire num
        // 2 * Math.tan(Math.PI / 6.) = 1.1547005383792515
        double y = (double) wire * 1.1547005383792515;
        if (layer % 2 == 1) {
            // Math.tan(Math.PI / 6.) = .5773502691896257
            y -= brickwallSign * 0.5773502691896257;
        }
        return y;
    }

    /**
     * Writes and returns a string containing the hit's important information.
     * @return hit's information encoded in a string
     */
    public String getInfo() {
        return "DC Hit:" +
               "\n    ID         : " + this.get_Id() +
               "\n    Sector     : " + this.get_Sector() +
               "\n    Superlayer : " + this.get_Superlayer() +
               "\n    Layer      : " + this.get_Layer() +
               "\n    Wire       : " + this.get_Wire() +
               "\n    TDC        : " + this.get_TDC();
    }

    /**
     * Writes and returns a string containing all of the hit's data.
     * @return hit's data encoded in a string
     */
    public String getDetailedInfo() {
        return "DC Hit " + this.get_Id() + ":" +
               "\n    Sector     : " + this.get_Sector() +
               "\n    Superlayer : " + this.get_Superlayer() +
               "\n    Layer      : " + this.get_Layer() +
               "\n    Wire       : " + this.get_Wire() +
               "\n    TDC        : " + this.get_TDC() +
               "\n    Cell size  : " + this.get_CellSize() +
               "\n    Doca err   : " + this.get_DocaErr() +
               "\n----------------------\n";
    }
}
