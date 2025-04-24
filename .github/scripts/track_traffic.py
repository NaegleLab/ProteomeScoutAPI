import requests
import json
import os
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt

def get_traffic_data(token):
    """Fetch and aggregate GitHub repository traffic data"""
    # Use environment variables set by GitHub Actions
    owner = os.environ.get('REPO_OWNER')
    repo = os.environ.get('REPO_NAME')
    
    headers = {
        'Authorization': f'token {token}',
        'Accept': 'application/vnd.github.v3+json'
    }
    
    base_url = f'https://api.github.com/repos/{owner}/{repo}'
    
    # Get current traffic data (14 days)
    views_response = requests.get(f'{base_url}/traffic/views', headers=headers)
    clones_response = requests.get(f'{base_url}/traffic/clones', headers=headers)
    referrers_response = requests.get(f'{base_url}/traffic/popular/referrers', headers=headers)
    
    # Process the data
    views = views_response.json()
    clones = clones_response.json()
    referrers = referrers_response.json() if referrers_response.status_code == 200 else []
    
    current_data = {
        'timestamp': datetime.now().isoformat(),
        'views': views,
        'clones': clones,
        'referrers': referrers
    }
    
    # Create directories if they don't exist
    os.makedirs('github_stats', exist_ok=True)
    os.makedirs('github_stats/raw', exist_ok=True)
    
    # Save raw data with timestamp
    raw_filename = f'github_stats/raw/{owner}-{repo}-{datetime.now().strftime("%Y%m%d")}.json'
    with open(raw_filename, 'w') as f:
        json.dump(current_data, f, indent=2)
    
    # Load and update aggregated data
    agg_filename = f'github_stats/{owner}-{repo}-aggregate.json'
    
    # Initialize aggregated data structure if it doesn't exist
    if not os.path.exists(agg_filename):
        aggregated_data = {
            'repository': f'{owner}/{repo}',
            'first_collected': datetime.now().isoformat(),
            'last_updated': datetime.now().isoformat(),
            'total_views': 0,
            'total_unique_views': 0,
            'total_clones': 0,
            'total_unique_clones': 0,
            'daily_views': {},
            'daily_clones': {},
            'referrers': {}
        }
    else:
        with open(agg_filename, 'r') as f:
            aggregated_data = json.load(f)
        aggregated_data['last_updated'] = datetime.now().isoformat()
    
    # Process views data
    for day_data in views.get('views', []):
        date = day_data.get('timestamp')[:10]  # Extract YYYY-MM-DD
        count = day_data.get('count', 0)
        uniques = day_data.get('uniques', 0)
        
        # Update daily views
        if date in aggregated_data['daily_views']:
            # Keep the highest count seen for this day
            existing_count = aggregated_data['daily_views'][date].get('count', 0)
            existing_uniques = aggregated_data['daily_views'][date].get('uniques', 0)
            
            aggregated_data['daily_views'][date] = {
                'count': max(count, existing_count),
                'uniques': max(uniques, existing_uniques)
            }
        else:
            # New date
            aggregated_data['daily_views'][date] = {
                'count': count,
                'uniques': uniques
            }
    
    # Process clones data
    for day_data in clones.get('clones', []):
        date = day_data.get('timestamp')[:10]  # Extract YYYY-MM-DD
        count = day_data.get('count', 0)
        uniques = day_data.get('uniques', 0)
        
        # Update daily clones
        if date in aggregated_data['daily_clones']:
            # Keep the highest count seen for this day
            existing_count = aggregated_data['daily_clones'][date].get('count', 0)
            existing_uniques = aggregated_data['daily_clones'][date].get('uniques', 0)
            
            aggregated_data['daily_clones'][date] = {
                'count': max(count, existing_count),
                'uniques': max(uniques, existing_uniques)
            }
        else:
            # New date
            aggregated_data['daily_clones'][date] = {
                'count': count,
                'uniques': uniques
            }
    
    # Process referrers
    for referrer in referrers:
        name = referrer.get('referrer', 'unknown')
        count = referrer.get('count', 0)
        uniques = referrer.get('uniques', 0)
        
        if name in aggregated_data['referrers']:
            aggregated_data['referrers'][name]['count'] += count
            aggregated_data['referrers'][name]['uniques'] += uniques
        else:
            aggregated_data['referrers'][name] = {
                'count': count,
                'uniques': uniques
            }
    
    # Recalculate totals
    total_views = sum(day['count'] for day in aggregated_data['daily_views'].values())
    total_unique_views = sum(day['uniques'] for day in aggregated_data['daily_views'].values())
    total_clones = sum(day['count'] for day in aggregated_data['daily_clones'].values())
    total_unique_clones = sum(day['uniques'] for day in aggregated_data['daily_clones'].values())
    
    aggregated_data['total_views'] = total_views
    aggregated_data['total_unique_views'] = total_unique_views 
    aggregated_data['total_clones'] = total_clones
    aggregated_data['total_unique_clones'] = total_unique_clones
    
    # Save updated aggregated data
    with open(agg_filename, 'w') as f:
        json.dump(aggregated_data, f, indent=2)
    
    # Print summary (will appear in GitHub Actions logs)
    print(f"\nTraffic Summary for {owner}/{repo}:")
    print(f"Current Views (14 days): {views.get('count', 0)} total, {views.get('uniques', 0)} unique")
    print(f"Current Clones (14 days): {clones.get('count', 0)} total, {clones.get('uniques', 0)} unique")
    print(f"All-time Views: {total_views} total")
    print(f"All-time Clones: {total_clones} total")
    print(f"Data collected since: {aggregated_data['first_collected']}")
    
    # Generate visualizations
    try:
        generate_traffic_report(aggregated_data, owner, repo)
    except Exception as e:
        print(f"Error generating visualizations: {e}")
    
    return aggregated_data

def generate_traffic_report(data, owner, repo):
    """Generate traffic reports and visualizations from aggregated data"""
    # Convert daily data to DataFrames for easier plotting
    views_data = []
    for date, stats in data['daily_views'].items():
        views_data.append({
            'date': date,
            'count': stats['count'],
            'uniques': stats['uniques']
        })
    
    clones_data = []
    for date, stats in data['daily_clones'].items():
        clones_data.append({
            'date': date,
            'count': stats['count'],
            'uniques': stats['uniques']
        })
    
    views_df = pd.DataFrame(views_data)
    clones_df = pd.DataFrame(clones_data)
    
    if not views_data or not clones_data:
        print("Not enough data to generate visualizations")
        return
    
    # Sort by date
    views_df['date'] = pd.to_datetime(views_df['date'])
    clones_df['date'] = pd.to_datetime(clones_df['date'])
    views_df = views_df.sort_values('date')
    clones_df = clones_df.sort_values('date')
    
    # Create visualizations
    plt.figure(figsize=(12, 8))
    
    # Views plot
    plt.subplot(2, 1, 1)
    plt.plot(views_df['date'], views_df['count'], 'b-', label='Total Views')
    plt.plot(views_df['date'], views_df['uniques'], 'g--', label='Unique Visitors')
    plt.title(f'GitHub Traffic for {owner}/{repo}')
    plt.ylabel('Views')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Clones plot
    plt.subplot(2, 1, 2)
    plt.plot(clones_df['date'], clones_df['count'], 'r-', label='Total Clones')
    plt.plot(clones_df['date'], clones_df['uniques'], 'm--', label='Unique Cloners')
    plt.xlabel('Date')
    plt.ylabel('Clones')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f'github_stats/{owner}-{repo}-traffic.png')
    
    # Create summary README file
    with open('github_stats/README.md', 'w') as f:
        f.write(f'# Traffic Statistics for {owner}/{repo}\n\n')
        f.write(f'*Last updated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}*\n\n')
        f.write(f'### Traffic Summary\n')
        f.write(f'- **Total Views**: {data["total_views"]}\n')
        f.write(f'- **Total Clones**: {data["total_clones"]}\n')
        f.write(f'- **Data Collection Started**: {data["first_collected"][:10]}\n\n')
        f.write(f'![Traffic Graph](/{owner}-{repo}-traffic.png)\n\n')
    
    # Create referrers visualization if there are referrers
    if data['referrers']:
        with open('github_stats/README.md', 'a') as f:
            f.write(f'### Top Referrers\n')
            f.write('| Referrer | Visits |\n')
            f.write('|----------|--------|\n')
            
            # Sort referrers by count
            sorted_referrers = sorted(
                [(name, stats['count']) for name, stats in data['referrers'].items()],
                key=lambda x: x[1],
                reverse=True
            )
            
            for name, count in sorted_referrers[:10]:  # Top 10 referrers
                f.write(f'| {name} | {count} |\n')

if __name__ == "__main__":
    # Get GitHub token from environment variables
    token = os.environ.get('GITHUB_TOKEN')
    if not token:
        raise ValueError("GitHub token not found in environment variables")
    
    get_traffic_data(token)