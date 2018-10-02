package org.jlab.service.kf;

import java.nio.ByteBuffer;

import org.jlab.io.base.DataEvent;
import org.jlab.rec.dc.track.Track;

public interface KFTrackCandDataEvent extends DataEvent {

    Track[] getTrack(String path);
    void    setTrack(String path, Track[] arr);
    void    appendTrack(String path, Track[] arr);
}
