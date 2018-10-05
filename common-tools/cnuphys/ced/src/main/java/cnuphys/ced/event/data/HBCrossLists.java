package cnuphys.ced.event.data;

public class HBCrossLists extends CrossLists {

	private static HBCrossLists _instance;

	private HBCrossLists() {
		super("HitBasedTrkg::HBCrossLists");
	}

	/**
	 * Public access to the singleton
	 * @return the singleton
	 */
	public static HBCrossLists getInstance() {
		if (_instance == null) {
			_instance = new HBCrossLists();
		}
		return _instance;
	}
}
