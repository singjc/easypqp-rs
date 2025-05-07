use assert_cmd::Command;

#[test]
fn test_cli_runs() {
    let mut cmd = Command::cargo_bin("easypqp").unwrap();
    cmd.arg("../test-data/config.json")
       .assert()
       .success();
}
