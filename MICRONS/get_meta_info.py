from caveclient import CAVEclient
from pathlib import Path 

client = CAVEclient()
## THIS WILL OPEN A WINDOW IN YOUR BROWSER
client.auth.setup_token(make_new=False)

## THIS NEEDS TO BE CHANGED TO THE TOKEN CREATED BY THE USER
client.auth.save_token(token="YOUR_TOKEN", overwrite=True)

# Access information
client = CAVEclient('minnie65_public')

# Desired save location, change to desired path.
save_dir = Path.home()

# Fetch orientation tuning information
orientation_info = client.materialize.query_table('digital_twin_properties_bcm_coreg_v4')
# Fetch neuron locations.
## Use the pt_position column for x-y-z coordinates.
location_info = client.materialize.query_table('coregistration_auto_phase3_fwd_apl_vess_combined_v2')

# Save data frames
location_info.to_csv(save_dir / "location.csv")
orientation_info.to_csv(save_dir / "orientation.csv")