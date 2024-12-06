import subprocess
import sys
import time
import os

def main():
    # ANSI color codes for the log
    BOLD = "\033[1m"
    RESET = "\033[0m"
    CYAN = "\033[1;36m"
    GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    RED = "\033[1;31m"
    MAGENTA = "\033[1;35m"

    CHECK_EMOJI = f"{GREEN}✔{RESET}"  # Green check mark for success
    ERROR_EMOJI = f"{RED}❌{RESET}"  # Red cross mark for errors
    module_dir = os.path.join(os.path.dirname(__file__), 'modules')

    # List of job scripts with their respective start and completion messages
    job_scripts = [
        {
            "script": os.path.join(module_dir, "install_packages.py"),
            "start_msg": "Installing required Python and R packages to run the analysis",
            "end_msg": "All Python and R packages are installed successfully.",
            "mandatory": True  # This step is mandatory to continue
        },
        {
            "script": os.path.join(module_dir, "data_preprocessing.py"),
            "start_msg": "Data Pre-Processing",
            "end_msg": "Data Pre-Processing Completed",
            "mandatory": True  # This step is mandatory after installing packages
        },
        {
            "script": os.path.join(module_dir, "best_recommended_k.py"),
            "start_msg": "Analyzing best K as per your choice of data and algorithm",
            "end_msg": "Analysis Completed",
            "mandatory": False  # Not mandatory, but will run if data preprocessing succeeds
        },
        {
            "script": os.path.join(module_dir, "clustering.py"),
            "start_msg": "Applying clustering algorithm to cluster the data",
            "end_msg": "Clustering Completed",
            "mandatory": True  # Mandatory for further analysis
        }
    ]

    # Additional analysis jobs with corresponding start and end messages
    analysis_jobs = {
        "1": {
            "script": os.path.join(module_dir, "survival_analysis.py"),
            "start_msg": "Performing Survival Analysis",
            "end_msg": "Survival Analysis Completed",
            "mandatory": False  # Example: Add 'True' for mandatory jobs
        },
        "2": {
            "script": os.path.join(module_dir, "mutation_analysis.py"),
            "start_msg": "Performing Mutation Analysis",
            "end_msg": "Mutation Analysis Completed",
            "mandatory": False  # Not mandatory
        },
        "3": {
            "script": os.path.join(module_dir, "stage_analysis.py"),
            "start_msg": "Performing Stage Analysis",
            "end_msg": "Stage Analysis Completed",
            "mandatory": False  # Not mandatory
        },
        "4": {
            "script": os.path.join(module_dir, "immune_analysis.py"),
            "start_msg": "Performing Immune Analysis",
            "end_msg": "Immune Analysis Completed",
            "mandatory": False  # Not mandatory
        }
    }

    def print_colored_analysis_complete():
        print(f"✅ {GREEN}Analysis Completed Successfully.{RESET}")

    def print_fancy_step(step, message):
        top_border = "✦" * 80
        bottom_border = "✦" * 80
        print(f"{top_border}")
        print(f"🌟 {MAGENTA}Step {step}:{RESET} {message}")
        print(bottom_border)

    def run_job(job, step):
        """Runs a Python script as a subprocess and logs its execution."""
        if step != 1:
            print("✨" * 20 + "✦" + "✨" * 20)
        print_fancy_step(step, job['start_msg'])

        try:
            # Log before running the job
            start_time = time.time()

            # Run the job script
            subprocess.check_call([sys.executable, job['script']])

            # Log after job completion
            end_time = time.time()
            execution_time = end_time - start_time
            print(f" {CHECK_EMOJI} {GREEN}{job['end_msg']}{RESET}")
            print(f"⌛ Execution Time: {execution_time:.2f} seconds")
            return True

        except subprocess.CalledProcessError as e:
            print(f"{RED}ERROR:{RESET} Error occurred in {job['start_msg']}: {e}")
            # print(f"Error occurred in {job['start_msg']}: {e}")
            return False
        except Exception as e:
            print(f"{RED}ERROR:{RED} An unexpected error occurred in {job['start_msg']}: {e}")
            return False

    def get_user_choice():
        """Displays a menu and returns the user's choice."""
        print("Select the further analysis you want to perform:")
        print("1. Survival Analysis")
        print("2. Mutation Analysis")
        print("3. Stage Analysis")
        print("4. Immune Analysis")
        print("5. Run All")
        return input("Enter the number of your choice: ")

    # Run initial jobs up to clustering
    step = 1
    install_success = run_job(job_scripts[0], step)
    if not install_success:
        print(f"{RED}FAILED: Packages Installation failed. Exiting...{RESET}")
        sys.exit(1)

    step += 1
    preprocessing_success = run_job(job_scripts[1], step)
    if not preprocessing_success:
        print(f"{RED}FAILED: Data Pre-Processing failed. Exiting...{RESET}")
        sys.exit(1)

    step += 1
    best_k_success = run_job(job_scripts[2], step)

    # Proceed with clustering analysis (must run even if Best K failed)
    step += 1
    clustering_success = run_job(job_scripts[3], step)
    if not clustering_success:
        print(f"{RED}FAILED: Clustering failed. Exiting...{RESET}")
        sys.exit(1)

    # Now ask user for further analysis steps
    choice = get_user_choice()

    # Run chosen analysis jobs based on user choice
    step += 1
    if choice == "5":
        # Run all further analysis jobs
        for key, analysis_job in analysis_jobs.items():
            if run_job(analysis_job, step):
                step += 1
    elif choice in analysis_jobs:
        # Run the selected analysis job
        run_job(analysis_jobs[choice], step)

    print_colored_analysis_complete()

# Call the main function
if __name__ == "__main__":
    main()
